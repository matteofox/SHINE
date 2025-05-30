#!/usr/bin/env python
# coding: utf-8
# AUTHORS: MF, DT
# VERSION: 1.1

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp

#import astropy modules
from astropy.io import fits
from astropy.wcs import WCS, utils
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.convolution import Gaussian1DKernel, Gaussian2DKernel, CustomKernel, interpolate_replace_nans
from astropy.table import Table, Column
from astropy.utils.exceptions import AstropyWarning

from scipy.ndimage import maximum_filter, median_filter
from astropy.convolution import convolve, convolve_fft
from astropy.stats import sigma_clipped_stats

from pathlib import Path

import argparse, textwrap

import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=AstropyWarning)

try:
    import cc3d
except:
    print('ERROR: This code requires the cc3d package. Please install the package using "pip install connected-components-3d" and try again')
    exit()

class SHINEWarning(UserWarning):
      pass
    
#Author MF
def filter_cube(cube, spatsmooth=2, specsig=0, isvar=False, usefftconv=False):
    
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
            for i in np.arange(cubsize[0]):
                cube[i, ...] = interpolate_replace_nans(cube[i, ...], tmpkern)
        else:
            label = 'data'
            normalize = True
            nan_treatment = 'interpolate'

        print('... Filtering the {} using XY-axis gaussian kernel of size {} pix'.format(label, spatsmooth))

        for i in np.arange(cubsize[0]):
            if usefftconv:
                SMcube[i, ...] = convolve_fft(cube[i, ...], spatkern, normalize_kernel=normalize,  nan_treatment=nan_treatment, allow_huge=True)
            else:
                SMcube[i, ...] = convolve(cube[i, ...], spatkern, normalize_kernel=normalize, nan_treatment=nan_treatment)
            
    if specsig > 0. and naxis==3:
        print('... Filtering the cube using Z-axis gaussian kernel of size {}'.format(specsig))
        speckern = Gaussian1DKernel(specsig)
        for i in np.arange(cubsize[1]):
          for j in np.arange(cubsize[2]):
            if usefftconv:
                SMcube[:,i,j] = convolve_fft(SMcube[:,i,j], speckern, normalize_kernel=normalize,  nan_treatment=nan_treatment, allow_huge=True)
            else:
                SMcube[:,i,j] = convolve(SMcube[:,i,j], speckern, normalize_kernel=normalize, nan_treatment=nan_treatment)
           

    elif specsig > 0. and naxis<3:   
        print('... Z-axis filtering requested on non-3D data. No spectral filtering will occurr.')
        
    if naxis==2:
       return SMcube[0,...]
    else:   
       return SMcube


def find_nan_edges(cube, extend=None):
    # If there are values identically zero, set them to NaN
    cube[cube == 0] = np.nan

    # Create a mask where all elements along the first axis are NaN
    mask = np.all(np.isnan(cube), axis=0)

    # Extend the mask if requested
    if extend is not None:
        if extend > 0:
            mask = maximum_filter(mask, size=extend)

    return mask
    
    

def threshold_cube(cube, var, threshold=0.5, maskedge=0, edge_ima=None):
    
    naxis = len(np.shape(cube))
    if naxis==2:
       cube = cube[np.newaxis,:]
       var = var[np.newaxis,:]
       print(f'... Thresholding in S/N the image with a threshold of {threshold}')
    else:
       print(f'... Thresholding in S/N the cube with a threshold of {threshold}')

    # Compute signal-to-noise ratio (S/N)
    snrcube = cube / np.sqrt(var)
    snrcube[np.logical_not(np.isfinite(snrcube))] = 0

    # Apply edge masking if requested
    if maskedge > 0:
        snrcube[:, :maskedge, :] = 0
        snrcube[:, :, :maskedge] = 0
        snrcube[:, -maskedge:, :] = 0
        snrcube[:, :, -maskedge:] = 0

    # Apply additional edge mask if provided
    if edge_ima is not None:
        edge_ima3d = edge_ima[np.newaxis, :]
        snrcube = snrcube * edge_ima3d

    # Create a boolean mask based on the threshold
    snrcube_bool = 1 * ((snrcube > threshold) & (np.isfinite(snrcube)))

    if naxis==2:
       return snrcube_bool[0,...], snrcube[0,...]
    else:
       return snrcube_bool, snrcube



def masking(data, mask):
    print('... Masking the threshold data')

    
    naxis = len(np.shape(data))
    if naxis==2:
         data = data[np.newaxis, :]
        
    # Adjust mask dimensions if necessary
    if np.shape(mask) == np.shape(data)[1:3]:
        mask3d = mask[np.newaxis, :]
    else:
        mask3d = mask

    # Identify bad regions where mask is greater than 0
    bad = mask3d[0, ...] > 0
    nz = np.shape(data)[0]

    # Replace bad regions with NaN
    for i in np.arange(nz):
        data[i, bad] = 0
        
    if naxis==2:
       return data[0,...]
    else:   
       return data


def masking_nan(data, mask):
    print('... Masking the data with nans')
    
    naxis = len(np.shape(data))
    if naxis==2:
       data = data[np.newaxis, :]

    # Adjust mask dimensions if necessary
    if np.shape(mask) == np.shape(data)[1:3]:
        mask3d = mask[np.newaxis, :]
    else:
        mask3d = mask
    
    # Identify bad regions where mask is greater than 0
    bad = mask3d[0, ...] > 0
    nz = np.shape(data)[0]

    # Replace bad regions with NaN
    for i in np.arange(nz):
        data[i, bad] = np.nan

    if naxis==2:
       return data[0,...]
    else:   
       return data
    
      


def subcube(cube=None, datahead=None, filename=None, pathcube=None, extcube=0, outdir='./', zmin=None, zmax=None, lmin=None, lmax=None, writesubcube=False, addname=''):
    
     
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
        subcube = cube[zmin:zmax, :, :]
    
    # If lmin and lmax are provided, calculate zmin and zmax
    elif lmin and lmax is not None:
        
        if lmin < min(wave):
            lmin = min(wave)
        if lmax > max(wave):
            lmax = max(wave)
            
        zmin = int(np.searchsorted(wave, lmin, side='right')-1)
        zmax = int(np.searchsorted(wave, lmax, side='right'))
        
        subcube = cube[zmin:zmax, :, :]
        
    else: 
        raise ValueError("Error: Please provide zmin/zmax or lmin/lmax.")

    #----------------------- SAVE THE OUTPUT ------------------------------- 
    print(f'... Selecting the cube between {zmin} and {zmax}')
    
    #update the header
    newhead['NAXIS3'] = zmax - zmin
    newhead['CRVAL3'] = wave[zmin]
    newhead['HISTORY'] = f'Cube selected between layer {zmin}-{zmax} using SHINE\'s subcube routine'
    
    if writesubcube:
        hduout = fits.PrimaryHDU(subcube, header = newhead)
        hduout.writeto(outdir+f'{filename}.SUBCUBE{addname}.fits', overwrite=True)
    #------------------------------------------------------------------------

        
    return subcube, newhead


#Author DT
def extract(cube, connectivity=26):

    print('... Extraction')

    #only 4,8 (2D) and 26, 18, and 6 (3D) are allowed
    
    naxis = len(np.shape(cube))
    if naxis == 2:
       if connectivity not in [4,8]:
          raise ValueError('Connectivity should be 4 or 8 for 2D images. Found {}.'.format(connectivity))
    elif naxis == 3:      
      if connectivity not in [6,18,26]:
          raise ValueError('Connectivity should be 6, 18, or 26 for 3D cubes. Found {}.'.format(connectivity))
          
    label_cube = cc3d.connected_components(cube, connectivity=connectivity)
    
    return label_cube


def generate_catalogue(labels, naxis=3):
    
    stats = cc3d.statistics(labels)
    
    IDs   = np.unique(labels)
    Nobj  = len(IDs)
    Nvox  = stats['voxel_counts']
    
    #Offset here is due to reference position for counting, the first pixel position should be 0.5 (center of first pixel)
    Xcent = stats['centroids'][:,-1]+0.5 
    Ycent = stats['centroids'][:,-2]+0.5
    #Evaluate Zcent only if 3D data
    if naxis==3:
       Zcent = stats['centroids'][:,-3]+0.5
    
    Xmin = []
    Xmax = []
    Ymin = []
    Ymax = []
    Zmin = []
    Zmax = []
    
    BBarray = stats['bounding_boxes']
        
    for ii in np.arange(Nobj):
        thisbb = BBarray[ii]
        
        Xmin.append(thisbb[-1].start)
        Xmax.append(thisbb[-1].stop) 
        Ymin.append(thisbb[-2].start)
        Ymax.append(thisbb[-2].stop)
        #Evaluate Zbb only if 3D data
        try:
         Zmin.append(thisbb[-3].start)
         Zmax.append(thisbb[-3].stop)
        except:
         pass
	
    Nspat = []
    
    Xmin = np.array(Xmin)
    Xmax = np.array(Xmax)
    Ymin = np.array(Ymin)
    Ymax = np.array(Ymax)
    Zmin = np.array(Zmin)
    Zmax = np.array(Zmax)        
        
    if naxis==3:
     for ind, tid in enumerate(IDs):
           pstamp = np.copy(labels[Zmin[ind]:Zmax[ind],Ymin[ind]:Ymax[ind],Xmin[ind]:Xmax[ind]])
            
           Nspat.append(calc_xyarea(pstamp,tid))
    
    Nspat = np.array(Nspat)
    Nz = Zmax-Zmin    #This should be right without a +1
    
    if naxis==2:
       table = Table([IDs, Nvox, Xcent, Ycent, Xmin, Xmax, Ymin, Ymax], \
            names=('ID', 'Npix', 'Xcent', 'Ycent', 'Xmin', 'Xmax', 'Ymin', 'Ymax'))
    elif naxis==3:       
       table = Table([IDs, Nvox, Nspat, Nz, Xcent, Ycent, Zcent, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax], \
            names=('ID', 'Nvox', 'Nspat', 'Nz', 'Xcent', 'Ycent', 'Zcent', 'Xmin', 'Xmax', 'Ymin', 'Ymax', 'Zmin', 'Zmax'))
    
    print('... Extraction DONE. {} objects found.'.format(len(IDs)))
    
    return table

def calc_xyarea(pstamp, tid):
    
    naxis = len(np.shape(pstamp))
    
    pstamp[pstamp!= tid] = 0
    if naxis>2:
       img = np.nanmax(pstamp, axis=0)
    else:
       img = pstamp
    
    return np.sum((img==tid))
    

def cleaning(catalogue, labels_out, mindz=1, maxdz=200, minvox=1, minarea=1):
    
    #If error is raised, it means dz is not available, i.e. the data is 2D
    try:
      dz = catalogue['Nz'] 
      keep = (catalogue['Nvox']>=minvox) & (dz>=mindz) & (dz<=maxdz) & (catalogue['Nspat']>=minarea)
      naxis = 3
    except:
      keep = (catalogue['Npix']>=np.nanmax([minvox, minarea]))
      naxis = 2
    
    remove = np.logical_not(keep)
    removeids = catalogue['ID'][remove]
    
    print('... Cleaning {} objects given the supplied criteria.'.format(len(removeids)))
    
    if naxis==2:
      for ind, tid in enumerate(removeids): 
        pstamp = np.copy(labels_out[catalogue['Ymin'][tid]:catalogue['Ymax'][tid],catalogue['Xmin'][tid]:catalogue['Xmax'][tid]])
        pstamp[pstamp==tid] = 0
        labels_out[catalogue['Ymin'][tid]:catalogue['Ymax'][tid],catalogue['Xmin'][tid]:catalogue['Xmax'][tid]] = pstamp
    elif naxis ==3:
      for ind, tid in enumerate(removeids): 
        pstamp = np.copy(labels_out[catalogue['Zmin'][tid]:catalogue['Zmax'][tid],catalogue['Ymin'][tid]:catalogue['Ymax'][tid],catalogue['Xmin'][tid]:catalogue['Xmax'][tid]])
        pstamp[pstamp==tid] = 0
        labels_out[catalogue['Zmin'][tid]:catalogue['Zmax'][tid],catalogue['Ymin'][tid]:catalogue['Ymax'][tid],catalogue['Xmin'][tid]:catalogue['Xmax'][tid]] = pstamp
    
    return catalogue[keep], labels_out
    

def add_wcs_struct(catalogue, datahead):
    
    wcs = WCS(datahead)
    naxis = datahead['NAXIS']
    
    skycoords_geom = utils.pixel_to_skycoord(catalogue['Xcent'], catalogue['Ycent'], wcs=wcs)
    skycoords_flxw = utils.pixel_to_skycoord(catalogue['XcentFL'], catalogue['YcentFL'], wcs=wcs)

    RA_column  = Column(data=np.round(skycoords_geom.ra.degree,6), name='RA_deg')   #in deg
    Dec_column = Column(data=np.round(skycoords_geom.dec.degree,6), name='DEC_deg') #in deg 

    catalogue.add_column(RA_column)
    catalogue.add_column(Dec_column)
        
    #Now build wave solution
    if naxis==3:
     try:
        dlam = datahead['CD3_3']
     except:
        dlam = datahead['CDELT3'] 
     finally:
        wave = datahead['CRVAL3']+ np.arange(datahead['NAXIS3'])*dlam
        
        Lambda_column = Column(data=np.interp(catalogue['Zcent'], np.arange(len(wave)), wave), name='Lambda') 
        catalogue.add_column(Lambda_column)   
  
    #Now move to flux weighted quantities
    RAFL_column  = Column(data=np.round(skycoords_flxw.ra.degree,6), name='RAFL_deg')   #in deg
    DecFL_column = Column(data=np.round(skycoords_flxw.dec.degree,6), name='DECFL_deg') #in deg 

    catalogue.add_column(RAFL_column)
    catalogue.add_column(DecFL_column)

    if naxis==3:
        
        LambdaFL_column = Column(data=np.interp(catalogue['ZcentFL'], np.arange(len(wave)), wave), name='LambdaFL') 
        catalogue.add_column(LambdaFL_column)   
     
    return catalogue


def compute_photometry(catalogue, cube, var, labelsCube):

    #As of April 2025, also compute flux weighted centroids in this function,
    #as for geometrical centroids, the first pixel position should be 0.5 (center of first pixel)
    
    naxis = len(np.shape(cube))
    
    flux     = np.zeros_like(catalogue['ID'], dtype=float)
    flux_err = np.zeros_like(catalogue['ID'], dtype=float)
    
    flwxcent = np.array(catalogue['Xmin'], dtype=float)+0.5
    flwycent = np.array(catalogue['Ymin'], dtype=float)+0.5
    if naxis==3:
       flwzcent = np.array(catalogue['Zmin'], dtype=float)+0.5
    
    for ind, tid in enumerate(catalogue['ID']): 
        
        if naxis==2:
          pstamplabel = labelsCube[catalogue['Ymin'][ind]:catalogue['Ymax'][ind],catalogue['Xmin'][ind]:catalogue['Xmax'][ind]]
          pstampcube  =       cube[catalogue['Ymin'][ind]:catalogue['Ymax'][ind],catalogue['Xmin'][ind]:catalogue['Xmax'][ind]]
          pstampvar   =        var[catalogue['Ymin'][ind]:catalogue['Ymax'][ind],catalogue['Xmin'][ind]:catalogue['Xmax'][ind]]
          
          okind = (pstamplabel==tid)
          valid_y, valid_x = np.nonzero(okind)
          
          
        elif naxis==3:
          pstamplabel = labelsCube[catalogue['Zmin'][ind]:catalogue['Zmax'][ind],catalogue['Ymin'][ind]:catalogue['Ymax'][ind],catalogue['Xmin'][ind]:catalogue['Xmax'][ind]]
          pstampcube  =       cube[catalogue['Zmin'][ind]:catalogue['Zmax'][ind],catalogue['Ymin'][ind]:catalogue['Ymax'][ind],catalogue['Xmin'][ind]:catalogue['Xmax'][ind]]
          pstampvar   =        var[catalogue['Zmin'][ind]:catalogue['Zmax'][ind],catalogue['Ymin'][ind]:catalogue['Ymax'][ind],catalogue['Xmin'][ind]:catalogue['Xmax'][ind]]

          okind = (pstamplabel==tid)
          valid_z, valid_y, valid_x = np.nonzero(okind)
         
        okind = (pstamplabel==tid)
        
        flux[ind] = np.nansum(pstampcube[okind])
        flux_err[ind] = np.sqrt(np.nansum(pstampvar[okind]))
        
        flwxcent[ind] += np.nansum(valid_x * pstampcube[okind])/flux[ind]
        flwycent[ind] += np.nansum(valid_y * pstampcube[okind])/flux[ind]
        if naxis==3:
           flwzcent[ind] += np.nansum(valid_z * pstampcube[okind])/flux[ind]
                
    Flux_column  = Column(data=flux, name='Flux')
    Flerr_column = Column(data=flux_err, name='Flux_err')
    SNR_column   = Column(data=flux/flux_err, name='Flux_SNR')
            
    catalogue.add_column(Flux_column)
    catalogue.add_column(Flerr_column)
    catalogue.add_column(SNR_column)
    
    XcentFL_column  = Column(data=flwxcent, name='XcentFL')
    YcentFL_column  = Column(data=flwycent, name='YcentFL')
    
    catalogue.add_column(XcentFL_column)
    catalogue.add_column(YcentFL_column)
 
    if naxis==3:
       ZcentFL_column  = Column(data=flwzcent, name='ZcentFL')
       catalogue.add_column(ZcentFL_column)
    
    return catalogue     



def compute_var(data):
    
    mean, median, std = sigma_clipped_stats(data, sigma = 3)
    vardata = (np.full_like(data, std))**2

    return vardata
    
    
 
def clean_clube(data, filtsize=7, rebinfac=40):
    
    nz, ny, nx = np.shape(data)

    data = np.ma.array(data, mask=np.isnan(data))
    
    zrebin = int(np.ceil(nz/rebinfac))
    
    contcube = np.zeros((zrebin, ny, nx))

    print(f'... Rebinning the cube')
    for ii in np.arange(zrebin):
        print(f'Slice: {ii}')
        zmin = rebinfac*ii
        zmax = min(nz, rebinfac*(ii+1))
        mean, median, std = sigma_clipped_stats(data[zmin:zmax,:,:], sigma = 3, axis=0, maxiters=3)   
        contcube[ii] = median

    print(f'... Filtering the rebinned cube')
    filtcube = median_filter(contcube, size=filtsize, axes=0)

    print(f'... Subtract the filtered cube')
    for ii in np.arange(zrebin):    
       
        zmin = rebinfac*ii
        zmax = min(nz, rebinfac*(ii+1))
        data[zmin:zmax,:,:] -= filtcube[ii]
    
    return data     



def runextraction(data, vardata, mask2d=None, mask2dpost=None, fmask3D=None, extdata=0, extvardata=0, \
                  snthreshold=2, maskspedge=0, spatsmooth=2, specsig=0, usefftconv=False, dovarsmooth=True, \
                  connectivity=26, mindz=1, maxdz=200, minvox = 1, minarea=1, \
                  zmin=None, zmax=None, lmin=None, lmax=None, outdir='./',  \
                  writelabels=False, writesmdata=False, writesmvar=False, writesmsnr=False, writesubcube=False, writevardata=False):

    
    hducube      = fits.open(data)
    hduhead      = hducube[extdata].header
    cubefilename = Path(data).stem

    # read or compute vardata
    try:
        vardata = float(vardata)
        varisnum = True
    except:
        varisnum = False
     
    if varisnum:  
        
        if vardata > 0:
            var = np.full_like(hducube[extdata].data, float(vardata))
            print('... Use the constant user-provided variance')
        elif vardata == -1:
            print('... Compute the variance with a sigma-clip algorithm')
            var = compute_var(hducube[extdata].data)
        else:
            raise ValueError('Error: vardata is not a valid number')

        varfilename = cubefilename+'_variance'
        
    else:
        hduvar = fits.open(vardata)
        var  = hduvar[extvardata].data
        varfilename  = Path(vardata).stem
        
    naxis = len(hducube[extdata].data.shape)
    
    #find edges
    edge_ima = find_nan_edges(hducube[extdata].data, extend=maskspedge)
       
          
    if mask2d is not None:
        mask2D = fits.open(mask2d)[0].data
        cube = masking_nan(hducube[extdata].data, mask2D)       
    else:
        cube = hducube[extdata].data
        
    
    #******************************************************************************************
    #step 0: create a subcube (in wavelength) of the original cube, if necessary and if data is 3D
    #******************************************************************************************
    if naxis==3:
      print('... Perform the extraction on a 3D data-cube')
      
      if zmin and zmax is not None:
        cube, newhduhead  = subcube(cube=cube, datahead=hduhead, filename=cubefilename, outdir=outdir, zmin=zmin, zmax=zmax, writesubcube=writesubcube)
        var, _   = subcube(cube=var,  datahead=hduhead, filename=varfilename,  outdir=outdir, zmin=zmin, zmax=zmax, writesubcube=writesubcube)
      elif lmin and lmax is not None:
        cube, newhduhead = subcube(cube=cube, datahead=hduhead, filename=cubefilename, outdir=outdir, lmin=lmin, lmax=lmax, writesubcube=writesubcube)
        var, _  = subcube(cube=var, datahead=hduhead, filename=varfilename, outdir=outdir, lmin=lmin, lmax=lmax, writesubcube=writesubcube)
      else:
        newhduhead = hduhead
        print('... Use the entire cube (without any selection in z direction)')
    
    elif naxis==2:
      print('... Perform the extraction on a 2D data-image')
      #Not elegant but necessary to conform to the specifications for 3D data
      newhduhead = hduhead

      if zmin or zmax is not None:
        warnings.warn('Zmin and/or Zmax arguments have been specified on a 2D input dataset. These arguments will be disregarded in the run.', SHINEWarning)
      elif lmin or lmax is not None:
         warnings.warn('Lmin and/or Lmax arguments have been specified on a 2D input dataset. These arguments will be disregarded in the run.', SHINEWarning)

    else:
      raise ValueError('Error: the data dimension is not correct. It can be a 2D or 3D dataset.')
          
       
    #******************************************************************************************
    #step 1: filtering the cube and the associated variance using a Gaussian kernel
    #******************************************************************************************
    cubeF = filter_cube(cube, spatsmooth=spatsmooth, specsig=specsig, usefftconv=usefftconv)
    if dovarsmooth:
       varF  = filter_cube(var,  spatsmooth=spatsmooth, specsig=specsig, usefftconv=usefftconv, isvar=True)
    else:
       varF  = var
    
    hducube.close()

    if not varisnum:  
        hduvar.close()
    
    
    #******************************************************************************************
    #step 2: extract the threshold cube
    #******************************************************************************************
    thcube, snrcube = threshold_cube(cubeF, varF, threshold=snthreshold, maskedge=maskspedge, edge_ima=1-edge_ima)
        
        
    #*********************************************************************************************
    #step 3: masking the voxels (needs masking again because nans are interpolated when filtering)
    #*********************************************************************************************
    if mask2dpost is not None:
        mask2Dpost = fits.open(mask2dpost)[0].data
        cubethresh = masking(thcube, mask2Dpost)
    else:
        cubethresh = thcube

        
    #*********************************************************************************************
    #step 4: grouping the voxels/pixels with chosen connectivity
    #*********************************************************************************************
    labels_out = extract(cubethresh, connectivity=connectivity)
   

    #*********************************************************************************************
    #step 5: generate a first version of the catalogue
    #*********************************************************************************************
    catalogue = generate_catalogue(labels_out, naxis=naxis)   
    
    
    #*********************************************************************************************
    #step 6: clean the catalogue and the labels map
    #*********************************************************************************************
    catalogue, labels_cln = cleaning(catalogue, labels_out, mindz=mindz, maxdz=maxdz, minvox=minvox, \
                                         minarea=minarea)
   
    
    #*********************************************************************************************
    #Step 7: perform photometry
    #*********************************************************************************************
    catalogue = compute_photometry(catalogue, cubeF, varF, labels_cln)
    
    
    #*********************************************************************************************
    #Step 8: add WCS in catalogue
    #*********************************************************************************************
    catalogue = add_wcs_struct(catalogue, newhduhead)
    
    
    #*********************************************************************************************
    #Write catalogue
    #*********************************************************************************************
    catalogue.write(outdir+f'/{cubefilename}.CATALOGUE_out.fits', overwrite=True)

    
    #*********************************************************************************************
    #Step 9: dump to the disk the requested products (note that the catalogue is always written)
    #*********************************************************************************************
    #Prepare header
    headout = newhduhead
    
    
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
    
    headout['XSMOOTH'] = xsig
    headout['YSMOOTH'] = ysig
    if naxis==3:
      headout['ZSMOOTH'] = specsig
    
    headout['SNTHRES'] = snthreshold
    
    
    if writelabels:
        hduout = fits.PrimaryHDU(labels_cln, header = headout)
        hduout.writeto(outdir+f'/{cubefilename}.LABELS_out.fits', overwrite=True)
            
    if writesmdata:
        hduout = fits.PrimaryHDU(cubeF, header = headout)
        hduout.writeto(outdir+f'/{cubefilename}.FILTER_out.fits', overwrite=True)
    
    if writesmvar:
        hduout = fits.PrimaryHDU(varF, header = headout)
        hduout.writeto(outdir+f'/{varfilename}.FILTER_out.fits', overwrite=True)
    
    if writesmsnr:
        hduout = fits.PrimaryHDU(snrcube, header = headout)
        hduout.writeto(outdir+f'/{cubefilename}.FILTERSNR_out.fits', overwrite=True)

    if writevardata:

        if varisnum:
            hduout = fits.PrimaryHDU(var, header = headout)
            hduout.writeto(outdir+f'/{varfilename}.fits', overwrite=True)
        else:
            print('... The variance has been provided as fits file. Not saving the file as an output.')
    
        

def main():
   
    parser = argparse.ArgumentParser(
    
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
    Extraction code for 2D/3D Data
    --------------------------------
    Authors: Matteo Fossati, Davide Tornotti
    
    '''))
    
    grpinp = parser.add_argument_group('Input control arguments') 
    
    grpinp.add_argument('data',            help='Path of the input data cube/image. Expected to be in extension 0, unless extdata is defined.')
    grpinp.add_argument('vardata',         help='Path of the variance cube/image.   Expected to be in extension 0, unless extvardata is defined. Instead, if a positive number is provided a constant variance is assumed. If -1, compute the variance from the data with a sigma-clipping algorithm.')
    grpinp.add_argument('--mask2d',        default= None,  help='Path of an optional two dimensional mask to be applied along the wave axis.')
    grpinp.add_argument('--mask2dpost',    default= None,  help='Path of an optional two dimensional mask to be applied after the smoothing along the wave axis.')
    grpinp.add_argument('--mask3d',        default= None,  help='Path of an optional three dimesional mask. Valid only for 3D data. NOT IMPLEMENTED YET.')
    grpinp.add_argument('--extdata',       default= 0, type=int, help='Specifies the HDU index in the FITS file cube to use for the data cube extraction.')
    grpinp.add_argument('--extvardata',    default= 0, type=int, help='Specifies the HDU index in the FITS file variance to use for the cube extraction.')
    grpinp.add_argument('--zmin',          default= None, type=int, help='If selecting the cube and the variance: initial pixel in z direction (from 0). Only valid for 3D data.')
    grpinp.add_argument('--zmax',          default= None, type=int, help='If selecting the cube and the variance: final pixel in z direction (from 0). Only valid for 3D data.')
    grpinp.add_argument('--lmin',          default= None, type=float, help='If selecting the cube and the variance: initial wavelength in z direction (in Angstrom). Only valid for 3D data.')
    grpinp.add_argument('--lmax',          default= None, type=float, help='If selecting the cube and the variance: final wavelength in z direction (in Angstrom). Only valid for 3D data.')
    
    grpext = parser.add_argument_group('Extraction arguments')
    
    grpext.add_argument('--snthreshold',     default= 2.,     type=float, help='The SNR of voxels to be included in the extraction.')
    grpext.add_argument('--spatsmooth',   default= 0.,     type=float, help='Gaussian Sigma of the spatial convolution kernel applied in X and Y.')
    grpext.add_argument('--spatsmoothX',  default= None,   help='Gaussian Sigma of the spatial convolution kernel applied in X. If set, this has priority over spatsmooth.')
    grpext.add_argument('--spatsmoothY',  default= None,   help='Gaussian Sigma of the spatial convolution kernel applied in Y. If set, this has priority over spatsmooth.')
    grpext.add_argument('--specsmooth',   default= 0.,     type=float, help='Gaussian Sigma of the spectra convolution kernel applied in Z/Lambda.')   
    grpext.add_argument('--usefftconv',   default= False,  type=bool, help='If True, use fft for convolution rather than the direct algorithm.')   
    grpext.add_argument('--dovarsmooth',  default= True,   type=bool, help='If False, do not apply the smoothing on the vardata.')
    grpext.add_argument('--connectivity', default= 26,     type=int, help='Voxel connectivity scheme to be used. Only 4,8 (2D) and 26, 18, and 6 (3D) are allowed.')   
    grpext.add_argument('--maskspedge',   default= 0,   type=int, help='Determines how much (in pixels) to expand the mask around the edges of the cube/image.')   
   
    grpcln = parser.add_argument_group('Cleaning arguments')
    
    grpcln.add_argument('--minvox',    default= 1,        type=int, help='Minimum number of connected voxels (3D)/pixels (2D) for a source to be in the final catalogue. For 2D data this argument has priority over minarea')
    grpcln.add_argument('--mindz',     default= 1,        type=int, help='Minimum number of connected voxels in spectral direction for a source to be in the final catalogue. Only valid for 3D data.')
    grpcln.add_argument('--maxdz',     default= 200,      type=int, help='Maximum number of connected voxels in spectral direction for a source to be in the final catalogue. Only valid for 3D data.')
    grpcln.add_argument('--minarea',   default= 1,        type=int, help='Minimum number of connected projected spatial voxels (3D)/pixels (2D) for a source to be in the final catalogue.')
    
    grpout = parser.add_argument_group('Output control arguments')
    
    grpout.add_argument('--outdir',              default='./',        help='Output directory path.')
    grpout.add_argument('--writelabels',         action='store_true', help='If set, write labels cube/image.')
    grpout.add_argument('--writesmdata',         action='store_true', help='If set, write the smoothed cube/image.')
    grpout.add_argument('--writesmvar',          action='store_true', help='If set, write the smoothed variance.')
    grpout.add_argument('--writesmsnr',          action='store_true', help='If set, write the S/N smoothed cube/image.')
    grpout.add_argument('--writesubcube',        action='store_true', help='If set and used, write the subcubes (cube and variance). Only valid for 3D data.')
    grpout.add_argument('--writevardata',        action='store_true', help='If vardata is a user-provided constant value or -1 (auto-computed from sigma-clipping) write the variance.')
    #...and more to come
        
    #Parse arguments  
    args = parser.parse_args()   
    
    #Manipulate some arguments  
    spatsig = (args.spatsmooth, args.spatsmooth)
    if args.spatsmoothX is not None:
        spatsig[0] = args.spatsmoothX
    if args.spatsmoothY is not None:
        spatsig[1] = args.spatsmoothY
    
    runextraction(args.data, args.vardata, \
    mask2d = args.mask2d, mask2dpost = args.mask2dpost, fmask3D = args.mask3d, extdata=args.extdata, extvardata=args.extvardata, \
    zmin = args.zmin, zmax=args.zmax, lmin=args.lmin, lmax=args.lmax, spatsmooth = spatsig, specsig=args.specsmooth, \
    usefftconv=args.usefftconv, dovarsmooth=args.dovarsmooth, snthreshold=args.snthreshold, connectivity = args.connectivity, \
    maskspedge=args.maskspedge, mindz = args.mindz, maxdz=args.maxdz, minvox=args.minvox, minarea=args.minarea, \
    outdir=args.outdir, writelabels=args.writelabels, writesmdata=args.writesmdata, \
    writesmvar=args.writesmvar, writesmsnr=args.writesmsnr, writevardata=args.writevardata)

if __name__ == "__main__":
   
   main()

