#!/usr/bin/env python
# coding: utf-8
# AUTHORS: MF, DT
# VERSION: 1.0

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

from scipy.ndimage import maximum_filter
from astropy.convolution import convolve, convolve_fft

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

    
#Author MF
def filter_cube(cube, spatsig=2, specsig=0, isvar=False, usefft=False):
    try:
        dummy = len(spatsig)
    except:
        spatsig = [spatsig]

    if len(spatsig) == 1:
        ysig = spatsig[0]
        xsig = spatsig[0]
    elif len(spatsig) > 1 and len(spatsig) <= 2:
        ysig = spatsig[1]
        xsig = spatsig[0]
    else:
        print('Error: the array of spatial smoothing kernels can have only 1 or 2 elements')
        return 0

    # Make a copy
    SMcube = np.copy(cube)
    # Get cube sizes
    cubsize = np.shape(cube)

    # Choose convolution method
    if usefft:
        myconv = convolve_fft
    else:
        myconv = convolve

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

        print('... Filtering the {} using XY-axis gaussian kernel of size {} pix'.format(label, spatsig))

        for i in np.arange(cubsize[0]):
            SMcube[i, ...] = myconv(cube[i, ...], spatkern, normalize_kernel=normalize, nan_treatment=nan_treatment)

    if specsig > 0.:
        print('... Filtering the cube using Z-axis gaussian kernel of size {}'.format(specsig))
        specsigel = Gaussian1DKernel(specsig)
        # Will be implemented

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

    return snrcube_bool, snrcube



def masking(cube_labels, mask):
    print('... Masking the labels cube')

    # Adjust mask dimensions if necessary
    if np.shape(mask) == np.shape(cube_labels)[1:3]:
        mask3d = mask[np.newaxis, :]
    else:
        mask3d = mask

    # Identify bad regions where mask is greater than 0
    bad = mask3d[0, ...] > 0
    nz = np.shape(cube_labels)[0]

    # Replace bad regions with NaN
    for i in np.arange(nz):
        cube_labels[i, bad] = 0
        
    return cube_labels


def masking_nan(cube, mask):
    print('... Masking the cube with nans')

    # Identify bad regions where mask is greater than 0
    bad = mask[0, ...] > 0
    nz = np.shape(cube)[0]

    # Replace bad regions with NaN
    for i in np.arange(nz):
        cube[i, bad] = np.nan

    return cube



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
            raise ValueError("Error: please provide the pathcube or data, header and filename.")
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
        print('... Please provide zmin/zmax or lmin/lmax')
    #------------------------------------------------------------------------
    
    
    
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

    # only 4,8 (2D) and 26, 18, and 6 (3D) are allowed

    label_cube = cc3d.connected_components(cube, connectivity=connectivity)

    return label_cube



def generate_catalogue(labels):
    
    stats = cc3d.statistics(labels)
    
    IDs   = np.unique(labels)
    Nobj  = len(IDs)
    Nvox  = stats['voxel_counts']
    Xcent = stats['centroids'][:,-1]
    Ycent = stats['centroids'][:,-2]
    Zcent = stats['centroids'][:,-3]
    
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
        Zmin.append(thisbb[-3].start)
        Zmax.append(thisbb[-3].stop)

    Nspat = []
    
    Xmin = np.array(Xmin)
    Xmax = np.array(Xmax)
    Ymin = np.array(Ymin)
    Ymax = np.array(Ymax)
    Zmin = np.array(Zmin)
    Zmax = np.array(Zmax)
    
        
    for ind, tid in enumerate(IDs):
        pstamp = np.copy(labels[Zmin[ind]:Zmax[ind],Ymin[ind]:Ymax[ind],Xmin[ind]:Xmax[ind]])
        Nspat.append(calc_xyarea(pstamp,tid))
    
    Nspat = np.array(Nspat)
    Nz = Zmax-Zmin    #This should be right without a +1
    
    table = Table([IDs, Nvox, Nspat, Nz, Xcent, Ycent, Zcent, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax], \
            names=('ID', 'Nvox', 'Nspat', 'Nz', 'Xcent', 'Ycent', 'Zcent', 'Xmin', 'Xmax', 'Ymin', 'Ymax', 'Zmin', 'Zmax'))
    
    print('... Extraction DONE. {} objects found.'.format(len(IDs)))
    
    return table

def calc_xyarea(pstamp, tid):
    
    pstamp[pstamp!= tid] = 0
    img = np.nanmax(pstamp, axis=0)
    
    return np.sum((img==tid))
    

def cleaning(catalogue, labels_out, mindz=1, maxdz=200, minvox=1, minarea=3):
    
    dz = catalogue['Nz'] 
        
    keep = (catalogue['Nvox']>=minvox) & (dz>=mindz) & (dz<=maxdz) & (catalogue['Nspat']>=minarea)
    remove = np.logical_not(keep)
    
    removeids = catalogue['ID'][remove]
    
    print('... Cleaning {} objects given the supplied criteria.'.format(len(removeids)))
    
    for ind, tid in enumerate(removeids): 
        pstamp = np.copy(labels_out[catalogue['Zmin'][tid]:catalogue['Zmax'][tid],catalogue['Ymin'][tid]:catalogue['Ymax'][tid],catalogue['Xmin'][tid]:catalogue['Xmax'][tid]])
        pstamp[pstamp==tid] = 0
        labels_out[catalogue['Zmin'][tid]:catalogue['Zmax'][tid],catalogue['Ymin'][tid]:catalogue['Ymax'][tid],catalogue['Xmin'][tid]:catalogue['Xmax'][tid]] = pstamp
    
    return catalogue[keep], labels_out
    

def add_wcs_struct(catalogue, datahead):
    
    wcs = WCS(datahead)
    
    skycoords = utils.pixel_to_skycoord(catalogue['Xcent'], catalogue['Ycent'], wcs=wcs)

    RA  = np.round(skycoords.ra.degree,6)  #in deg
    Dec = np.round(skycoords.dec.degree,6) #in deg 
    
    RA_column  = Column(data=RA, name='RA_deg')
    Dec_column = Column(data=Dec, name='DEC_deg')
    
    #Now build wave solution
    try:
        try:
            dlam = datahead['CD3_3']
        except:
            dlam = datahead['CDELT3'] 
                     
        wave = datahead['CRVAL3']+ np.arange(datahead['NAXIS3'])*dlam
        
        Lambda = np.interp(catalogue['Zcent'], np.arange(len(wave)), wave)
    except:
        Lambda = np.ones_like(catalogue['Zcent'])
        
    
    Lambda_column = Column(data=Lambda, name='Lambda')  

    catalogue.add_column(RA_column)
    catalogue.add_column(Dec_column)
    catalogue.add_column(Lambda_column)  
    
    return catalogue


def compute_photometry(catalogue, cube, var, labelsCube):
    
    flux     = np.zeros_like(catalogue['ID'], dtype=float)
    flux_err = np.zeros_like(catalogue['ID'], dtype=float)
    
    for ind, tid in enumerate(catalogue['ID']): 
        
        pstamplabel = labelsCube[catalogue['Zmin'][ind]:catalogue['Zmax'][ind],catalogue['Ymin'][ind]:catalogue['Ymax'][ind],catalogue['Xmin'][ind]:catalogue['Xmax'][ind]]
        pstampcube  =       cube[catalogue['Zmin'][ind]:catalogue['Zmax'][ind],catalogue['Ymin'][ind]:catalogue['Ymax'][ind],catalogue['Xmin'][ind]:catalogue['Xmax'][ind]]
        pstampvar   =        var[catalogue['Zmin'][ind]:catalogue['Zmax'][ind],catalogue['Ymin'][ind]:catalogue['Ymax'][ind],catalogue['Xmin'][ind]:catalogue['Xmax'][ind]]
        
        okind = (pstamplabel==tid)
        
        flux[ind] = np.nansum(pstampcube[okind])
        flux_err[ind] = np.sqrt(np.nansum(pstampvar[okind]))
                
    Flux_column  = Column(data=flux, name='Flux')
    Flerr_column = Column(data=flux_err, name='Flux_err')
    SNR_column   = Column(data=flux/flux_err, name='Flux_SNR')
                
    catalogue.add_column(Flux_column)
    catalogue.add_column(Flerr_column)
    catalogue.add_column(SNR_column)
    
    return catalogue     


def runextraction(fcube, fvariance, fmask2D=None, fmask2Dpost=None, fmask3D=None, extcub=0, extvar=0, \
                  SNthreshold=2, maskspedge=0, spatsig=2, specsig=0, usefft=False, connectivity=26, \
                  mindz=1, maxdz=200, minvox = 1, minarea=3, zmin=None, zmax=None, lmin=None, lmax=None, outdir='./', \
                  writelabels=False, writesmcube=False, writesmvar=False, writesmsnrcube=False, writesubcube=False):

    
    hducube      = fits.open(fcube)
    hduvar       = fits.open(fvariance)
    hduhead      = hducube[extcub].header
    cubefilename = Path(fcube).stem
    varfilename  = Path(fvariance).stem
    
    #find edges
    edge_ima        = find_nan_edges(hducube[extcub].data, extend=maskspedge)
       
          
    if fmask2D is not None:
        mask2D = fits.open(fmask2D)[0].data
        cube = masking_nan(hducube[extcub].data,   mask2D)
        #var  = masking_nan(hduvar[extvar].data,    mask2D)
        var  = hduvar[extvar].data
    else:
        cube = hducube[extcub].data
        var  = hduvar[extvar].data        
        
    
    #******************************************************************************************
    #step 0: create a subcube (in wavelength) of the original cube, if necessary
    #******************************************************************************************
    if zmin and zmax is not None:
        cube, newhduhead  = subcube(cube=cube, datahead=hduhead, filename=cubefilename, outdir=outdir, zmin=zmin, zmax=zmax, writesubcube=writesubcube)
        var, _   = subcube(cube=var,  datahead=hduhead, filename=varfilename,  outdir=outdir, zmin=zmin, zmax=zmax, writesubcube=writesubcube)
    elif lmin and lmax is not None:
        cube, newhduhead = subcube(cube=cube, datahead=hduhead, filename=cubefilename, outdir=outdir, lmin=lmin, lmax=lmax, writesubcube=writesubcube)
        var, _  = subcube(cube=var, datahead=hduhead, filename=varfilename, outdir=outdir, lmin=lmin, lmax=lmax, writesubcube=writesubcube)
    else:
        newhduhead = hduhead
        print('... Use the entire cube (without any selection in z direction)')
    
    
    #******************************************************************************************
    #step 1: filtering the cube and the associated variance using a Gaussian kernel
    #******************************************************************************************
    cubeF = filter_cube(cube, spatsig=spatsig, specsig=specsig, usefft=usefft)
    varF  = filter_cube(var,  spatsig=spatsig, specsig=specsig, usefft=usefft, isvar=True)
    
    hducube.close()
    hduvar.close()
    
    
    #******************************************************************************************
    #step 2: extract the threshold cube
    #******************************************************************************************
    thcube, snrcube = threshold_cube(cubeF, varF, threshold=SNthreshold, maskedge=maskspedge, edge_ima=1-edge_ima)
        
        
    #*********************************************************************************************
    #step 3: masking the voxels (needs masking again because nans are interpolated when filtering)
    #*********************************************************************************************
    if fmask2Dpost is not None:
        mask2Dpost = fits.open(fmask2Dpost)[0].data
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
    catalogue = generate_catalogue(labels_out)   
    
    
    #*********************************************************************************************
    #step 6: clean the catalogue and the labels map
    #*********************************************************************************************
    catalogue, labels_cln = cleaning(catalogue, labels_out, mindz=mindz, maxdz=maxdz, minvox=minvox, \
                                         minarea=minarea)

    
    #*********************************************************************************************
    #Step 7: add WCS in catalogue
    #*********************************************************************************************
    catalogue = add_wcs_struct(catalogue, newhduhead)
    
    
    #*********************************************************************************************
    #Step 8: perform photometry
    #*********************************************************************************************
    catalogue = compute_photometry(catalogue, cubeF, varF, labels_cln)
    
    
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
        dummy = len(spatsig)
    except:
        spatsig = [spatsig]

    if len(spatsig) == 1:
        ysig = spatsig[0]
        xsig = spatsig[0]
    elif len(spatsig) > 1 and len(spatsig) <= 2:
        ysig = spatsig[1]
        xsig = spatsig[0]
    
    headout['XSMOOTH'] = xsig
    headout['YSMOOTH'] = ysig
       
    headout['ZSMOOTH'] = specsig
    headout['SNTHRES'] = SNthreshold
    
    
    if writelabels:
        hduout = fits.PrimaryHDU(labels_cln, header = headout)
        hduout.writeto(outdir+f'/{cubefilename}.LABELS_out.fits', overwrite=True)
            
    if writesmcube:
        hduout = fits.PrimaryHDU(cubeF, header = headout)
        hduout.writeto(outdir+f'/{cubefilename}.FILTER_out.fits', overwrite=True)
    
    if writesmvar:
        hduout = fits.PrimaryHDU(varF, header = headout)
        hduout.writeto(outdir+f'/{varfilename}.FILTER_out.fits', overwrite=True)
    
    if writesmsnrcube:
        hduout = fits.PrimaryHDU(snrcube, header = headout)
        hduout.writeto(outdir+f'/{cubefilename}.SNRcubeF_out.fits', overwrite=True)
        

def main():
   
    parser = argparse.ArgumentParser(
    
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
    Extraction code for 2D/3D Data
    --------------------------------
    Authors: Matteo Fossati, Davide Tornotti
    
    '''))
    
    grpinp = parser.add_argument_group('Input control arguments') 
    
    grpinp.add_argument('cube',        help='Path of the input datacube. Expected to be in extension 0, unless extcub is defined.')
    grpinp.add_argument('varcube',     help='Path of the variance cube. Expected to be in extension 0, unless extvar is defined.')
    grpinp.add_argument('--mask2d',    default= None,  help='Path of an optional two dimensional mask to be applied along the wave axis.')
    grpinp.add_argument('--mask2dpost',default= None,  help='Path of an optional two dimensional mask to be applied after the smoothing along the wave axis.')
    grpinp.add_argument('--mask3d',    default= None,  help='Path of an optional three dimesional mask. NOT IMPLEMENTED YET.')
    grpinp.add_argument('--extcub',    default= 0, type=int, help='Specifies the HDU index in the FITS file cube to use for the data cube extraction.')
    grpinp.add_argument('--extvar',    default= 0, type=int, help='Specifies the HDU index in the FITS file variance to use for the cube extraction.')
    grpinp.add_argument('--zmin',      default= None, type=int, help='If selecting the cube and the variance: initial pixel in z direction (from 0).')
    grpinp.add_argument('--zmax',      default= None, type=int, help='If selecting the cube and the variance: final pixel in z direction (from 0).')
    grpinp.add_argument('--lmin',      default= None, type=float, help='If selecting the cube and the variance: initial wavelength in z direction (in Angstrom).')
    grpinp.add_argument('--lmax',      default= None, type=float, help='If selecting the cube and the variance: final wavelength in z direction (in Angstrom).')
    
    
    grpext = parser.add_argument_group('Extraction arguments')
    
    grpext.add_argument('--snthresh',     default= 2.,     type=float, help='The SNR of voxels to be included in the extraction.')
    grpext.add_argument('--spatsmooth',   default= 0.,     type=float, help='Gaussian Sigma of the spatial convolution kernel applied in X and Y.')
    grpext.add_argument('--spatsmoothX',  default= None,   help='Gaussian Sigma of the spatial convolution kernel applied in X. If set, this has priority over spatsmooth.')
    grpext.add_argument('--spatsmoothY',  default= None,   help='Gaussian Sigma of the spatial convolution kernel applied in Y. If set, this has priority over spatsmooth.')
    grpext.add_argument('--specsmooth',   default= 0.,     type=float, help='Gaussian Sigma of the spectra convolution kernel applied in Lambda. NOT IMPLEMENTED YET.')   
    grpext.add_argument('--usefftconv',   default= False,  type=bool, help='If True, use fft for convolution rather than the direct algorithm.')   
    grpext.add_argument('--connectivity', default= 26,     type=int, help='Voxel connectivity scheme to be used. Only 4,8 (2D) and 26, 18, and 6 (3D) are allowed.')   
    grpext.add_argument('--maskspedge',   default= None,   type=int, help='Determines how much (in pixels) to expand the mask around the edges of the cube/image.')   
   
    grpcln = parser.add_argument_group('Cleaning arguments')
    
    grpcln.add_argument('--minvox',    default= 1,        type=int, help='Minimum number of connected voxels for a source to be in the final catalogue.')
    grpcln.add_argument('--mindz',     default= 1,        type=int, help='Minimum number of connected voxels in spectral direction for a source to be in the final catalogue.')
    grpcln.add_argument('--maxdz',     default= 200,      type=int, help='Maximum number of connected voxels in spectral direction for a source to be in the final catalogue.')
    grpcln.add_argument('--minarea',   default= 3,        type=int, help='Minimum number of connected voxels in projected spatial direction for a source to be in the final catalogue.')
    
    grpout = parser.add_argument_group('Output control arguments')
    
    grpout.add_argument('--outdir',              default='./',        help='Output directory path.')
    grpout.add_argument('--writelabels',         action='store_true', help='If set, write labels cube.')
    grpout.add_argument('--writesmcube',         action='store_true', help='If set, write the smoothed cube.')
    grpout.add_argument('--writesmvar',          action='store_true', help='If set, write the smoothed variance.')
    grpout.add_argument('--writesmsnrcube',      action='store_true', help='If set, write the S/N smoothed cube.')
    grpout.add_argument('--writesubcube',        action='store_true', help='If set and used, write the subcubes (cube and variance).')
    
    #...and more to come
        
    #Parse arguments  
    args = parser.parse_args()   
    
    #Manipulate some arguments  
    spatsig = (args.spatsmooth, args.spatsmooth)
    if args.spatsmoothX is not None:
        spatsig[0] = args.spatsmoothX
    if args.spatsmoothY is not None:
        spatsig[1] = args.spatsmoothY
    
    #if args.maskspedge==None:
    #    args.maskspedge=int(5*spatsig[0])
    
    runextraction(args.cube, args.varcube, \
    fmask2D = args.mask2d, fmask2Dpost = args.mask2dpost, fmask3D = args.mask3d, extcub=args.extcub, extvar=args.extvar, \
    spatsig = spatsig, specsig=args.specsmooth, usefft=args.usefftconv, SNthreshold=args.snthresh, \
    maskspedge=args.maskspedge, mindz = args.mindz, maxdz=args.maxdz, minvox=args.minvox, minarea=args.minarea, \
    outdir=args.outdir, writelabels=args.writelabels, writesmcube=args.writesmcube, \
    writesmvar=args.writesmvar, writesmsnrcube=args.writesmsnrcube)

if __name__ == "__main__":
   
   main()

