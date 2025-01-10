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
from pathlib import Path
import os

warnings.filterwarnings("ignore", category=RuntimeWarning)




def Make_Im_SHINE(cube, labelsCube=None, Id=[-1], extcub=0, extlabels=0, itype='mean', outdir='./', writeout=False, addname='', nsl=-2, nsladd=0):
    
    
    #------------------ READ THE CUBE AND SET THE DATA ---------------------
    hducube  = fits.open(cube)
    hduhead  = hducube[extcub].header
    filename = Path(cube).stem

    try: 
        deltal  = hduhead['CD3_3']
    except:
        deltal  = hduhead['CDELT3']

    try:
      pixsize = np.abs(np.sqrt(hduhead['CD1_1']**2+hduhead['CD1_2']**2))*3600
    except:  
      pixsize = np.abs(hduhead['CDELT1'])*3600 #assuming the same for X and Y (transform to arcsec)
        
    #read the data
    cube   = np.nan_to_num(hducube[extcub].data)
    #-------------------------------------------------------------------- 
    
    

    #-------------- USE ID INFORMATION IF PROVIDED ----------------------
    #if the ID cube is used produce image with the sources
    if labelsCube is not None:
        hdulabel = fits.open(labelsCube)
        labels = hdulabel[extlabels].data

        # If all objects requested make a full list
        if Id[0] == -1:
            Id = np.unique(labels)
            Id = Id[1:]
         
        
        
        # if only one object use this ID
        if len(Id) == 1:
            
            keepvox = (labels == Id[0])

            if (nsl!=-2) and (nsl >= -1):
                # Ensure nsl works only with one single ID
                z_indices = np.where(keepvox.any(axis=(1, 2)))[0]  # Get all z indices for the given ID

                if len(z_indices) == 0:
                    raise ValueError(f"No voxels found for ID {Id[0]} in the label cube.")
                
                if nsl==-1:
                    central_layer = int(np.round(np.mean(z_indices)))  # Find the central layer
                else:
                    central_layer = nsl


                # Define the range of layers to include (central +/- nsladd)
                min_layer = max(0, central_layer - nsladd)
                max_layer = min(cube.shape[0], central_layer + nsladd + 1)

                # Create a mask to keep the voxels in the range of layers
                layer_mask = np.zeros(cube.shape, dtype=bool)
                layer_mask[min_layer:max_layer, :, :] = True


                # Combine the layer mask with the ID mask
                combined_mask = keepvox | layer_mask

                # Zero out voxels not in the combined mask
                cube[~combined_mask] = 0
                
                
            elif nsl==-2:
 
                # If nsl is -2, just keep the voxels for the given ID
                trimvox = np.logical_not(keepvox)
                cube[trimvox] = 0         
                
            else:
                raise ValueError("Please, provide a valid nsl (default = -2)")
                 
             
        else:
            if (nsl!=-2) and (nsl >= -1):
                raise ValueError("If nsl != -2, you must use a single ID to create the image.")

            # Combine all selected IDs if more than one
            keepvox = (labels == Id[0])
            for i in range(1, len(Id)):
                keepvox |= (labels == Id[i])

            trimvox = np.logical_not(keepvox)
            cube[trimvox] = 0
                       
    #----------------------------------------------------------------------- 
        
        
        
    #-------------- PRODUCE THE IMAGE GIVEN THE METHOD ---------------------
    if itype == 'flux':
        image = np.nansum(cube, axis=0)
        image    = (image*deltal/pixsize**2)*0.01 #to have units in 10^-18 erg s^-1 cm^-2 arcsec^-1
     
    elif itype =='mean':
        image = np.nanmean(cube, axis=0)
    
    elif itype =='median':
        image = np.nanmedian(cube, axis=0)
        
    else:
        print('Error: provide a valid image type (flux, mean, median)')
        return
             
    image[image==0]=np.nan
    #----------------------------------------------------------------------- 
    
    
    if writeout:
        
        hduout_img = fits.PrimaryHDU(image[np.newaxis,:,:])
        
        #Edit the header
        headout = hduout_img.header
        
        copykeys = ['OBJECT', 'CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2', 'CUNIT1', 'CUNIT2', 'CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2']
        for key in copykeys:
            if key in hduhead:  
                headout[key] = hduhead[key]
            else:
                if key == 'CDELT1':  
                    if 'CD1_D1' in hduhead:
                        headout[key] = hduhead['CD1_D1']
                elif key == 'CDELT2':  
                    if 'CD2_D2' in hduhead:
                        headout[key] = hduhead['CD2_D2']
        
        if itype == 'flux':
            headout['BUNIT'] = '1e-18 erg cm^-2 s^-1 arcsec^-2'
        else:
            headout['BUNIT'] = '1e-20 erg cm^-2 s^-1 Angstrom^-1'
        
        headout['HISTORY'] = f'2D image using Make_Im_SHINE with method {itype}'
        
        #Paste the header and save
        hduout_img.header = headout
        hduout_img.writeto(outdir+f'{filename}.IMAGE{addname}.fits', overwrite=True)
            
    return image



def main():
   
    parser = argparse.ArgumentParser(
    
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
    Make surface brightness images of 3D data using the output of SHINE or create narrow band images. 
    --------------------------------
    Authors:  Davide Tornotti, Matteo Fossati
    
    '''))
    
    grpinp = parser.add_argument_group('Input control arguments') 
    
    grpinp.add_argument('cube',            help='Path of the input datacube. Expected to be in extension 0, unless extcub is defined.')
    grpinp.add_argument('labelsCube',      default = None, help='Path of the cube with labels. Expected to be in extension 0, unless extlabels is defined.')
    grpinp.add_argument('--Id',            nargs="+",  default=[-1], type=int, help='The Ids of the grouped voxels to be used for the surface brightness image extraction. If a list is passed, all the valid IDs will be stacked into the final image. If [-1] (Default) stacks all the Ids in the labelsCube.')
    grpinp.add_argument('--itype',         default= 'flux', help='Type of image to produce: <flux>, <mean> or <median>')
    grpinp.add_argument('--extcub',        default= 0, type=int, help='Extension including science data')
    grpinp.add_argument('--extlabels',     default= 0, type=int, help='Extension including labels data')
    grpinp.add_argument('--nsl',           default=-2, type=int,  help='Noise layers: select this layer associated with the object in the image. If -1 it selects the mean layer. Use this only if len(Id)=1. Default is -2 -> no noise layers.')
    grpinp.add_argument('--nsladd',        default= 0, type=int, help='Noise layers: how many layers collapse adjacent to the selected one.')
    
    
    
    grpout = parser.add_argument_group('Output control arguments')
    
    grpout.add_argument('--outdir',        default='./',        help='Output directory path')
    grpout.add_argument('--writeout',      action='store_true', help='If set, write flux image and associated variance image')
    grpinp.add_argument('--addname',       default='', help='Optional suffix to append to the base name of the output file. ' 
                         'If set, this string will be added after the default filename (e.g., "flux" or other predefined parts), '
                         'before the file extension. Useful for distinguishing different versions or types of output files.')
    
    #Parse arguments  
    args = parser.parse_args()   
    
    if args.Id[0] >=0:
        print('The following IDs will be collapsed in the image: {}'.format(args.Id))
    else:
        print('All the extracted labels will be collapsed in the image')
  
      
    Make_Im_SHINE(args.cube, args.labelsCube, args.Id, itype = args.itype, extcub=args.extcub, extlabels=args.extlabels, outdir=args.outdir, writeout = args.writeout, addname = args.addname, nsl=args.nsl, nsladd=args.nsladd)
    
    
if __name__ == "__main__":
    main()
