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
from astropy.convolution import convolve, convolve_fft, Gaussian1DKernel, Gaussian2DKernel, CustomKernel
from astropy.table import Table, Column

import argparse, textwrap
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)




def Make_Im_SHINE(cube, variance, labelsCube, Id, extcub=0, extvar=0, extlabels=0, outdir='./', writeout=False):
    
    hducube  = fits.open(cube)
    hduvar   = fits.open(variance)
    hdulabel = fits.open(labelsCube)
    
    hduhead  = hducube[extcub].header
    
    #read the data
    cube   = np.nan_to_num(hducube[extcub].data)
    var    = np.nan_to_num(hduvar[extvar].data)
    labels = hdulabel[extlabels].data
    
    #If all objects requested make a full list
    if Id[0] == -1:
       Id = np.unique(labels)
       Id = Id[1:]
    
    lamdain = hduhead['CRVAL3']
    try: 
      deltal  = hduhead['CD3_3']
    except:
      deltal  = hduhead['CDELT3']
      
    pixsize = np.abs(hduhead['CDELT2'])*3600 #assuming the same for CDELT1 (transform to arcsec)
    
    keepvox = (labels==Id[0])
    if len(Id)>1:
       for i in range(1, len(Id)):
          keepvox = keepvox | (labels==Id[i])
    
    trimvox = np.logical_not(keepvox)
    cube[trimvox] = 0
    var[trimvox] = 0
    
    image = np.nansum(cube, axis=0)
    varimage = np.nansum(var, axis=0)         
    
    image    = (image*deltal/pixsize**2)*0.01 #to have units in 10^-18 erg s^-1 cm^-2 arcsec^-1
    varimage = (varimage*deltal/pixsize**2)*0.01 
    
    image[image==0]=np.nan
    varimage[varimage==0]=np.nan
    
    
    if writeout:
        
        hduout_img = fits.PrimaryHDU(image)
        hduout_var = fits.PrimaryHDU(varimage)
        
        #Edit the header
        headout = hduout_img.header
        
        copykeys = ['OBJECT','WCSAXES','CRPIX1','CRPIX2','CDELT1','CDELT2','CUNIT1','CUNIT2','CTYPE1','CTYPE2','CRVAL1','CRVAL2','LONPOLE','LATPOLE','MJDREF','RADESYS']
        for key in copykeys:
            headout[key] = hduhead[key]
        
        headout['BUNIT'] = '1e-18 erg cm^-2 s^-1 arcsec^-2'
        headout['HISTORY'] = '2D image from SHINE extraction'
        
        #Paste the header
        hduout_img.header = headout
        hduout_var.header = headout
        
        hduout_img.writeto(outdir+'C2image.fits', overwrite=True)
        hduout_var.writeto(outdir+'C2varimage.fits', overwrite=True)
        
        
    return image, varimage


if __name__ == "__main__":
   
   
    parser = argparse.ArgumentParser(
    
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
    Make surface brightness images of 3D data using the output of SHINE 
    Units of output image 1e-18 erg s^-1 cm^-2 arcsec^-1
    --------------------------------
    Authors:  Davide Tornotti, Matteo Fossati
    
    '''))
    
    grpinp = parser.add_argument_group('Input control arguments') 
    
    grpinp.add_argument('cube',            help='Path of the input datacube. Expected to be in extension 0, unless extcub is defined.')
    grpinp.add_argument('varcube',         help='Path of the variance cube. Expected to be in extension 0, unless extvar is defined.')
    grpinp.add_argument('labelsCube',      help='Path of the variance cube. Expected to be in extension 0, unless extvar is defined.')
    grpinp.add_argument('--Id',            nargs="+",  default=[-1], type=int, help='The Ids of the grouped voxels to be used for the surface brightness image extraction. If a list is passed, all the valid IDs will be stacked into the final image')
    grpinp.add_argument('--extcub',        default= 0, type=int, help='Extension including science data')
    grpinp.add_argument('--extvar',        default= 0, type=int, help='Extension including variance data')
    grpinp.add_argument('--extlabels',     default= 0, type=int, help='Extension including labels data')
    
    grpout = parser.add_argument_group('Output control arguments')
    
    grpout.add_argument('--outdir',        default='./',        help='Output directory path')
    grpout.add_argument('--writeout',      action='store_true', help='If set, write flux image and associated variance image')
    
    #Parse arguments  
    args = parser.parse_args()   
    
    if args.Id[0] >=0:
       print('The following IDs will be collapsed in the image: {}'.format(args.Id))
    else:
       print('All the extracted labels will be collapsed in the image')
  
      
    Make_Im_SHINE(args.cube, args.varcube, args.labelsCube, args.Id, extcub=args.extcub, extvar=args.extvar, extlabels=args.extlabels, outdir=args.outdir, writeout = args.writeout)
    
    
