#!/usr/bin/env python
# coding: utf-8
# AUTHORS: MF, DT
# VERSION: 2.0
#
# shine_utils.py
# --------------
# General-purpose utility functions for the SHINE pipeline.
#
# Functions specific to the emitter-extraction workflow (covariance
# estimation, catalogue building, etc.) have been moved to the
# ``shine.Find_Em_SHINE`` sub-package.


import warnings

import numpy as np
from scipy.ndimage import median_filter

from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.io import fits
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel, CustomKernel, interpolate_replace_nans
import concurrent.futures
from pathlib import Path

warnings.filterwarnings("ignore", category=RuntimeWarning)


# =============================================================================
# 1.  clean_clube — continuum subtraction
# =============================================================================

def clean_clube(data, filtsize=7, rebinfac=40):
    """Subtract a smooth continuum from a 3-D spectroscopic cube.

    The cube is first re-binned spectrally by collapsing groups of
    ``rebinfac`` layers into a single median-stack slice. The resulting
    low-resolution continuum cube is then smoothed along the spectral axis
    with a sliding-window median filter of width ``filtsize``. The smoothed
    continuum is finally subtracted layer-by-layer from the original cube.

    Parameters
    ----------
    data : numpy.ndarray, shape (nz, ny, nx)
        Input 3-D data cube.  NaN values are treated as masked pixels via
        :class:`numpy.ma.MaskedArray`.
    filtsize : int, optional
        Width (in re-binned spectral slices) of the median filter applied to
        the low-resolution continuum cube. Default is 7.
    rebinfac : int, optional
        Number of original spectral layers collapsed into each continuum
        slice. Default is 40.

    Returns
    -------
    data : numpy.ma.MaskedArray, shape (nz, ny, nx)
        Continuum-subtracted cube.  The mask mirrors the NaN pattern of the
        input.
    """
    nz, ny, nx = np.shape(data)

    data = np.ma.array(data, mask=np.isnan(data))

    zrebin = int(np.ceil(nz / rebinfac))

    contcube = np.zeros((zrebin, ny, nx))

    print(f'... Rebinning the cube into {zrebin} continuum slices '
          f'(rebinfac={rebinfac})')
          
    def _process_slice(ii):
        zmin = rebinfac * ii
        zmax = min(nz, rebinfac * (ii + 1))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            clipped = sigma_clip(data[zmin:zmax, :, :], sigma=3, maxiters=3, axis=0)
            filled = clipped.filled(np.nan)
            return ii, np.nanmedian(filled, axis=0)

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(_process_slice, ii) for ii in range(zrebin)]
        for future in concurrent.futures.as_completed(futures):
            ii, median = future.result()
            contcube[ii] = median

    print(f'... Filtering the continuum cube with a spectral median filter '
          f'(filtsize={filtsize})')
    filtcube = median_filter(contcube, size=filtsize, axes=0)

    print('... Subtracting the smooth continuum')
    for ii in np.arange(zrebin):
        zmin = rebinfac * ii
        zmax = min(nz, rebinfac * (ii + 1))
        data[zmin:zmax, :, :] -= filtcube[ii]

    return data


# =============================================================================
# 2.  filter_cube — spatial and spectral smoothing
# =============================================================================

def _process_interpolate_nans(args):
    i, data_slice, kern = args
    return i, interpolate_replace_nans(data_slice, kern)

def _process_convolve(args):
    i, data_slice, kern, normalize, nan_treatment, usefft = args
    if usefft:
        return i, convolve_fft(data_slice, kern, normalize_kernel=normalize, nan_treatment=nan_treatment, allow_huge=True)
    else:
        return i, convolve(data_slice, kern, normalize_kernel=normalize, nan_treatment=nan_treatment)

def filter_cube(cube, spatsmooth=2, specsig=0, isvar=False, usefftconv=False):
    """Apply a 2-D or 3-D Gaussian smoothing kernel to a data cube or variance.

    Parameters
    ----------
    cube : numpy.ndarray
        Input 2-D image or 3-D data cube.
    spatsmooth : float or list of float, optional
        Gaussian sigma (in pixels) for the spatial smoothing. If a single float
        is provided, the same sigma is applied in X and Y. If a list of two
        floats is provided, they correspond to [sig_x, sig_y]. Default is 2.
    specsig : float, optional
        Gaussian sigma (in pixels) for the spectral smoothing (Z axis). Only
        valid for 3-D data. Default is 0 (disabled).
    isvar : bool, optional
        If True, the smoothing is tailored for variance data (the kernel is
        squared and not normalized to unity). Default is False.
    usefftconv : bool, optional
        If True, use FFT-based convolution instead of direct convolution.
        Default is False.

    Returns
    -------
    SMcube : numpy.ndarray
        Smoothed data cube or image with the same shape as the input.
    """
    try:
        dummy = len(spatsmooth)
    except:
        spatsmooth = [spatsmooth]

    if len(spatsmooth) == 1:
        ysig = spatsmooth[0]
        xsig = spatsmooth[0]
    elif len(spatsmooth) > 1 and len(spatsmooth) <= 2:
        ysig = spatsmooth[1]
        xsig = spatsmooth[0]
    else:
        raise ValueError('Error: the spatial smoothing kernels can be an integer or an array with 1 or 2 elements')

    # Make a copy
    SMcube = np.copy(cube)
    # Get cube sizes
    cubsize = np.shape(cube)
    naxis = len(cubsize)
    
    if naxis==2:
       SMcube = SMcube[np.newaxis,:]
       cube   = cube[np.newaxis,:]
       cubsize = np.shape(cube)

    if ysig > 0. and xsig > 0.:

        if specsig == 0:

            # ----------------------------------------------------------------
            # this is the spatial 2D smoothing case (valid for 2D and 3D data)
            # ----------------------------------------------------------------
            spatkern = Gaussian2DKernel(xsig, ysig, x_size=int(6 * xsig + 1), y_size=int(6 * ysig + 1))

            if isvar:
                # Variance requires a special treatment because the kernel cannot be normalized to unity
                label = 'variance'
                nan_treatment = 'fill'
                normalize = False
    
                # Make a custom Kernel
                spatkern = CustomKernel((spatkern.array) ** 2)
    
                # Interpolate NaNs with ad-hoc kernel
                print('... Interpolating NaNs in Variance Data')
                tmpkern = Gaussian2DKernel(xsig, ysig, x_size=int(6 * xsig + 1), y_size=int(6 * ysig + 1))
                with concurrent.futures.ProcessPoolExecutor() as executor:
                    args_list = [(i, cube[i, ...], tmpkern) for i in np.arange(cubsize[0])]
                    futures = [executor.submit(_process_interpolate_nans, args) for args in args_list]
                    for future in concurrent.futures.as_completed(futures):
                        i, result = future.result()
                        cube[i, ...] = result
            else:
                label = 'data'
                normalize = True
                nan_treatment = 'interpolate'
    
            print('... Filtering the {} using XY-axis gaussian kernel of size xsig={}, ysig={} pix'.format(label, xsig, ysig))
    
            with concurrent.futures.ProcessPoolExecutor() as executor:
                args_list = [(i, cube[i, ...], spatkern, normalize, nan_treatment, usefftconv) for i in np.arange(cubsize[0])]
                futures = [executor.submit(_process_convolve, args) for args in args_list]
                for future in concurrent.futures.as_completed(futures):
                    i, result = future.result()
                    SMcube[i, ...] = result
                    
        elif specsig > 0. and naxis==3:

            # -----------------------------------------------------------------------
            # this is the spatial+spectral 3D smoothing case (valid only for 3D data)
            # -----------------------------------------------------------------------
            spatspeckern = Gaussian3DKernel(xsig, ysig, specsig, xsize=int(6 * xsig + 1), ysize=int(6 * ysig + 1), zsize=int(6 * specsig + 1) )
    
            if isvar:
                # Variance requires a special treatment because the kernel cannot be normalized to unity
                label = 'variance'
                nan_treatment = 'fill'
                normalize = False
    
                # Make a custom Kernel
                spatspeckern = CustomKernel((spatspeckern) ** 2)
    
                # Interpolate NaNs with ad-hoc kernel
                print('... Interpolating NaNs in Variance Data')
                tmpkern = Gaussian3DKernel(xsig, ysig, specsig, xsize=int(6 * xsig + 1), ysize=int(6 * ysig + 1), zsize=int(6 * specsig + 1))
                cube = interpolate_replace_nans(cube, tmpkern)
                
            else:
                label = 'data'
                normalize = True
                nan_treatment = 'interpolate'
    
            print('... Filtering the {} using XYZ-axis gaussian kernel of size xsig={}, ysig={} and zsig={} pix'.format(label, xsig, ysig, specsig))
    
            if usefftconv:
                SMcube = convolve_fft(cube, spatspeckern, normalize_kernel=normalize,  nan_treatment=nan_treatment, allow_huge=True)
            else:
                SMcube = convolve(cube, spatspeckern, normalize_kernel=normalize, nan_treatment=nan_treatment)
        else:   
            raise ValueError('... Z-axis filtering requested on non-3D data.')

    elif ysig == 0. and xsig == 0.:
        if naxis==2:
           return SMcube[0,...]
        else:   
           return SMcube
    else:
        raise ValueError('... Invalid xsig and ysig. They must be > 0')
            
    if naxis==2:
       return SMcube[0,...]
    else:   
       return SMcube


def Gaussian3D(xstd, ystd, zstd, xmean=0, ymean=0, zmean=0):
    amplitude = 1/( (2*np.pi)**(3/2)*xstd*ystd*zstd)
    def gaussian(x,y,z):
        f = amplitude*np.exp( -((x-xmean)**2/(2*xstd**2) + (y-ymean)**2/(2*ystd**2) + (z-zmean)**2/(2*zstd**2)))
        return f
    return gaussian


def Gaussian3DKernel(xstd, ystd, zstd, xsize, ysize, zsize):
    g = Gaussian3D(xstd, ystd, zstd)
    x = np.arange(- (xsize // 2), (xsize // 2)+1) 
    y = np.arange(- (ysize // 2), (ysize // 2)+1) 
    z = np.arange(- (zsize // 2), (zsize // 2)+1) 
    zz, yy, xx = np.meshgrid(z, y, x, indexing='ij')
    kernel_array = g(xx, yy, zz)
    kernel_array /= np.sum(kernel_array)
    return kernel_array


# =============================================================================
# 3.  subcube — spectral sub-cube extraction
# =============================================================================

def subcube(cube=None, datahead=None, filename=None, pathcube=None, extcube=0, outdir='./', zmin=None, zmax=None, lmin=None, lmax=None, writesubcube=False, addname=''):
    """Extract a spectral sub-cube by pixel index (Z axis) or wavelength range.

    Parameters
    ----------
    cube : numpy.ndarray, optional
        Input 3-D data cube array.
    datahead : astropy.io.fits.Header, optional
        FITS header of the input cube (used to extract wavelength WCS keywords).
    filename : str, optional
        Stem name of the file (used for writing the output file).
    pathcube : str, optional
        Path to the FITS file of the cube. If provided, `cube`, `datahead`, and
        `filename` are read directly from it.
    extcube : int, optional
        HDU extension index to read from `pathcube`. Default is 0.
    outdir : str, optional
        Output directory path. Default is './'.
    zmin, zmax : int, optional
        Spectral layer pixel index range (0-based) to extract.
    lmin, lmax : float, optional
        Wavelength range (in Å) to extract.
    writesubcube : bool, optional
        If True, write the extracted sub-cube to a FITS file. Default is False.
    addname : str, optional
        Suffix to append to the output filename.

    Returns
    -------
    subcube : numpy.ndarray
        The extracted 3-D sub-cube.
    newhead : astropy.io.fits.Header
        The updated FITS header reflecting the new spectral range.
    """
    #----------------------- PRELIMINARY CHECKS ----------------------------
    #if set the patcube read the data from here 
    if pathcube is not None:
        cube     = fits.open(pathcube)[extcube].data
        datahead = fits.open(pathcube)[extcube].header
        filename = Path(pathcube).stem
    else:
        #if the pathcube is not provided check the data, datahead, filename
        if cube is None or datahead is None or filename is None:
            raise ValueError("Error: please provide the pathcube or data+header+filename.")
    #-----------------------------------------------------------------------
    
    #----------------------- SELECT THE CUBE -------------------------------
    # assuming zmin starts from 0 and/or lmin from the minimum wavelength
    
    # Build wavelengths 
    try:
        dlam = datahead['CD3_3']
    except:
        dlam = datahead['CDELT3'] 
    
    wave = datahead['CRVAL3'] + np.arange(datahead['NAXIS3']) * dlam
    
    newhead = datahead.copy()
    
    # Extract the subcube based on zmin and zmax
    if zmin and zmax is not None:
        subcube_data = cube[zmin:zmax, :, :]
    
    # If lmin and lmax are provided, calculate zmin and zmax
    elif lmin and lmax is not None:
        
        if lmin < min(wave):
            lmin = min(wave)
        if lmax > max(wave):
            lmax = max(wave)
            
        zmin = int(np.searchsorted(wave, lmin, side='right')-1)
        zmax = int(np.searchsorted(wave, lmax, side='right'))
        
        subcube_data = cube[zmin:zmax, :, :]
        
    else: 
        raise ValueError("Error: Please provide zmin/zmax or lmin/lmax.")

    #----------------------- SAVE THE OUTPUT ------------------------------- 
    print(f'... Selecting the cube between {zmin} and {zmax}')
    
    #update the header
    newhead['NAXIS3'] = zmax - zmin
    newhead['CRVAL3'] = wave[zmin]
    newhead['HISTORY'] = f'Cube selected between layer {zmin}-{zmax} using SHINE\'s subcube routine'
    
    if writesubcube:
        hduout = fits.PrimaryHDU(subcube_data, header = newhead)
        hduout.writeto(outdir+f'{filename}.SUBCUBE{addname}.fits', overwrite=True)
    #------------------------------------------------------------------------
        
    return subcube_data, newhead

