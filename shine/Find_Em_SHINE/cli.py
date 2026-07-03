#!/usr/bin/env python
# coding: utf-8
# AUTHORS: MF, DT
# VERSION: 2.0
#
# Find_Em_SHINE/cli.py
# --------------------
# Command-line interface for the emitter extraction pipeline using config.ini.

import os
import sys
import argparse
import textwrap
import configparser
import numpy as np

# SHINE imports
from .extraction import extract
from .covariance import estimate_empirical_covariance as covariance
from .catalogue import build_emitter_catalogue as build_em_catalog


# =============================================================================
# Default Template Config
# =============================================================================
TEMPLATE_CONFIG = """; Find_Em_SHINE configuration file
; Generate a default template with: Find_Em_SHINE --init-config filename.ini

[FILES]
; Science data cube FITS file path (required)
fcube = /path/to/Datacube.fits
; Variance cube FITS file path (required)
fvar = /path/to/Varcube.fits
; Output directory for SHINE extraction products
outdir_extraction = ./extraction/
; Output directory for covariance estimation products
outdir_covariance = ./covariance/
; Output directory for the final catalogue and cutout images
outdir_catalogue = ./emitters/
; HDU extension index for the science data cube
extdata = 0
; HDU extension index for the variance cube
extvar = 0
; Path to an optional pre-smoothing 2-D mask (1 = bad)
mask2d = None
; Path to an optional post-smoothing 2-D mask (1 = bad)
mask2dpost = None

[CONTINUUM_SUBTRACTION]
; Enable continuum subtraction before extraction
do_continuum_sub = False
; Spectral re-binning factor (layers collapsed per continuum slice)
rebinfac = 40
; Median filter width (in re-binned slices) along spectral axis
filtsize = 7

[EXTRACTION]
; S/N threshold for voxel inclusion
snthreshold = 3.0
; Spatial Gaussian smoothing kernel sigma (pixels)
spatsmooth = 2.0
; Spectral Gaussian smoothing kernel sigma (pixels)
specsmooth = 0.0
; Voxel connectivity scheme (6, 18, or 26)
connectivity = 26
; Pixels to mask around the field edges
maskspedge = 0
; Minimum spectral extent (layers) per source
mindz = 3
; Maximum spectral extent (layers) per source
maxdz = 50
; Minimum total connected voxels per source
minvox = 27
; Minimum projected spatial area (pixels) per source
minarea = 9
; Layer index range selection (e.g. zmin=40, zmax=100) or None
zmin = None
zmax = None
; Wavelength range selection (in Angstrom) or None
lmin = None
lmax = None

[COVARIANCE]
; Number of spectral wavelength segments
nsegments = 500
; Aperture half-sizes (pixels) to probe (comma-separated list)
allsizes = 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30
; Spectral depth (layers) of each random aperture sample
dl = 4
; Valid samples to collect per aperture size per segment
nsamples = 10000
; Pixel scale in arcseconds per pixel (MUSE nominal is 0.2)
pixel_scale = 0.2
; Degree of polynomial fit to the covariance vs aperture size
fitdeg = 2
; Path to a 2-D source mask (or None to rely on automatic sep detection)
mask_source_cov = None
; Generate a 2-panel diagnostic PDF plot
plot = True
; Print detailed progress information
verbose = True

[CATALOGUE]
; Target redshift for velocity-offset computation
target_z = None
; Rest-frame wavelength of the emission line (e.g., Ly-alpha at 1216.0)
rest_line = None
; Maximum velocity offset (km/s) to keep a source in the catalog
vel_cut = None
; S/N thresholds defining confidence classes (comma-separated, highest first)
sncut = 7.0, 5.0
; Generate per-source cutout images (cutout FITS files)
checkimg = True
; Spatial padding (pixels) around each cutout postage-stamp
padding = 50
; Path to the median data cube (optional, for SNR_med checks)
fcube_median = None
; Paths to half-exposure data cubes for even/odd checks (optional)
fcube_odd = None
fcube_even = None
; Associated variance cubes for median and half-exposure cubes
fcube_median_var = None
fcube_odd_var = None
fcube_even_var = None
; Path to the unsmoothed data cube for 1-D spectrum extraction (requires mypython)
; Defaults to the input cube (fcube) if do_continuum_sub=True and this is None.
fcube_for_spectra = None
; Path to a 2-D continuum source mask for overlap flagging
fsource_img = None
; Path to a Marz-format redshift catalogue for continuum sources
marzred = None
; Max fractional even/odd S/N difference per class (comma-separated)
delta_eo_sncut = 0.5, 0.5
; Minimum even and odd S/N per class (comma-separated)
sn_eo_cut = 3.0, 3.0
; Minimum fraction of voxels within a 5x5x5 bounding box
fracnpix = None
"""


# =============================================================================
# Helper function for parsing types
# =============================================================================
def _parse_val(val, val_type=str):
    if val is None:
        return None
    val_str = str(val).strip()
    if val_str.lower() in ('none', 'null', ''):
        return None
    
    if val_type == bool:
        return val_str.lower() in ('true', 'yes', '1', 'on')
    
    if val_type == 'list_int':
        return [int(x.strip()) for x in val_str.split(',') if x.strip()]
    
    if val_type == 'tuple_float':
        return tuple(float(x.strip()) for x in val_str.split(',') if x.strip())
    
    return val_type(val_str)


# =============================================================================
# Command-line Execution
# =============================================================================
def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        Find_Em_SHINE — Emitter Extraction & Cataloguing CLI
        ----------------------------------------------------
        Authors: Davide Tornotti, Matteo Fossati
        
        Runs the full 3-step pipeline (extraction, covariance, catalogue)
        using parameters defined in a config.ini file.
        ''')
    )

    parser.add_argument(
        'config',
        nargs='?',
        default=None,
        help='Path to the config.ini file.'
    )
    
    parser.add_argument(
        '--init-config',
        dest='init_config',
        default=None,
        help='Initialize a template config.ini file at the specified path.'
    )

    args = parser.parse_args()

    # Handle --init-config flag
    if args.init_config is not None:
        filepath = args.init_config
        try:
            with open(filepath, 'w') as f:
                f.write(TEMPLATE_CONFIG)
            print(f"Created template configuration file: {filepath}")
            sys.exit(0)
        except Exception as e:
            print(f"ERROR: Could not write template configuration file: {e}")
            sys.exit(1)

    # Make sure we have a config path if not initializing
    if args.config is None:
        parser.print_help()
        sys.exit(0)

    config_path = args.config
    if not os.path.isfile(config_path):
        print(f"ERROR: Configuration file not found: {config_path}")
        sys.exit(1)

    print(f"Reading configuration from: {config_path}")
    
    # Parse the config file
    config = configparser.ConfigParser()
    config.read(config_path)

    # ------------------------------------------------------------------
    # Section: FILES
    # ------------------------------------------------------------------
    try:
        fcube = _parse_val(config.get('FILES', 'fcube'), str)
        fvar  = _parse_val(config.get('FILES', 'fvar'), str)
    except configparser.NoSectionError as e:
        print(f"ERROR: Missing section in config file: {e}")
        sys.exit(1)
    except configparser.NoOptionError as e:
        print(f"ERROR: Missing required option: {e}")
        sys.exit(1)

    if not fcube or not fvar:
        print("ERROR: fcube and fvar must be specified in [FILES] section.")
        sys.exit(1)

    outdir_extraction = _parse_val(config.get('FILES', 'outdir_extraction', fallback='./extraction/'), str)
    outdir_covariance = _parse_val(config.get('FILES', 'outdir_covariance', fallback='./covariance/'), str)
    outdir_catalogue  = _parse_val(config.get('FILES', 'outdir_catalogue', fallback='./emitters/'), str)

    extdata    = _parse_val(config.get('FILES', 'extdata', fallback='0'), int)
    extvar     = _parse_val(config.get('FILES', 'extvar', fallback='0'), int)
    mask2d     = _parse_val(config.get('FILES', 'mask2d', fallback=None), str)
    mask2dpost = _parse_val(config.get('FILES', 'mask2dpost', fallback=None), str)

    # ------------------------------------------------------------------
    # Section: CONTINUUM_SUBTRACTION
    # ------------------------------------------------------------------
    do_continuum_sub = _parse_val(config.get('CONTINUUM_SUBTRACTION', 'do_continuum_sub', fallback='False'), bool)
    rebinfac         = _parse_val(config.get('CONTINUUM_SUBTRACTION', 'rebinfac', fallback='40'), int)
    filtsize         = _parse_val(config.get('CONTINUUM_SUBTRACTION', 'filtsize', fallback='7'), int)

    # ------------------------------------------------------------------
    # Section: EXTRACTION
    # ------------------------------------------------------------------
    snthreshold  = _parse_val(config.get('EXTRACTION', 'snthreshold', fallback='2.0'), float)
    spatsmooth   = _parse_val(config.get('EXTRACTION', 'spatsmooth', fallback='2.0'), float)
    specsmooth   = _parse_val(config.get('EXTRACTION', 'specsmooth', fallback='0.0'), float)
    connectivity = _parse_val(config.get('EXTRACTION', 'connectivity', fallback='26'), int)
    maskspedge   = _parse_val(config.get('EXTRACTION', 'maskspedge', fallback='20'), int)
    mindz        = _parse_val(config.get('EXTRACTION', 'mindz', fallback='1'), int)
    maxdz        = _parse_val(config.get('EXTRACTION', 'maxdz', fallback='200'), int)
    minvox       = _parse_val(config.get('EXTRACTION', 'minvox', fallback='1'), int)
    minarea      = _parse_val(config.get('EXTRACTION', 'minarea', fallback='1'), int)
    zmin         = _parse_val(config.get('EXTRACTION', 'zmin', fallback=None), int)
    zmax         = _parse_val(config.get('EXTRACTION', 'zmax', fallback=None), int)
    lmin         = _parse_val(config.get('EXTRACTION', 'lmin', fallback=None), float)
    lmax         = _parse_val(config.get('EXTRACTION', 'lmax', fallback=None), float)

    # ------------------------------------------------------------------
    # Section: COVARIANCE
    # ------------------------------------------------------------------
    nsegments       = _parse_val(config.get('COVARIANCE', 'nsegments', fallback='500'), int)
    allsizes_str    = config.get('COVARIANCE', 'allsizes', fallback='2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30')
    allsizes        = _parse_val(allsizes_str, 'list_int')
    dl              = _parse_val(config.get('COVARIANCE', 'dl', fallback='4'), int)
    nsamples        = _parse_val(config.get('COVARIANCE', 'nsamples', fallback='10000'), int)
    pixel_scale     = _parse_val(config.get('COVARIANCE', 'pixel_scale', fallback='0.2'), float)
    fitdeg          = _parse_val(config.get('COVARIANCE', 'fitdeg', fallback='2'), int)
    mask_source_cov = _parse_val(config.get('COVARIANCE', 'mask_source_cov', fallback=None), str)
    plot_covariance = _parse_val(config.get('COVARIANCE', 'plot', fallback='True'), bool)
    verbose_cov     = _parse_val(config.get('COVARIANCE', 'verbose', fallback='True'), bool)

    # ------------------------------------------------------------------
    # Section: CATALOGUE
    # ------------------------------------------------------------------
    target_z         = _parse_val(config.get('CATALOGUE', 'target_z', fallback=None), float)
    rest_line        = _parse_val(config.get('CATALOGUE', 'rest_line', fallback=None), float)
    vel_cut          = _parse_val(config.get('CATALOGUE', 'vel_cut', fallback=None), float)
    sncut            = _parse_val(config.get('CATALOGUE', 'sncut', fallback='7.0,5.0'), 'tuple_float')
    checkimg         = _parse_val(config.get('CATALOGUE', 'checkimg', fallback='True'), bool)
    padding          = _parse_val(config.get('CATALOGUE', 'padding', fallback='50'), int)
    fcube_median     = _parse_val(config.get('CATALOGUE', 'fcube_median', fallback=None), str)
    fcube_odd        = _parse_val(config.get('CATALOGUE', 'fcube_odd', fallback=None), str)
    fcube_even       = _parse_val(config.get('CATALOGUE', 'fcube_even', fallback=None), str)
    fcube_median_var = _parse_val(config.get('CATALOGUE', 'fcube_median_var', fallback=None), str)
    fcube_odd_var    = _parse_val(config.get('CATALOGUE', 'fcube_odd_var', fallback=None), str)
    fcube_even_var   = _parse_val(config.get('CATALOGUE', 'fcube_even_var', fallback=None), str)
    fcube_for_spectra = _parse_val(config.get('CATALOGUE', 'fcube_for_spectra', fallback=None), str)
    if fcube_for_spectra is None:
        # Fallback to fcube_orig for backward compatibility
        fcube_for_spectra = _parse_val(config.get('CATALOGUE', 'fcube_orig', fallback=None), str)

    if fcube_for_spectra is None and do_continuum_sub:
        fcube_for_spectra = fcube

    fsource_img      = _parse_val(config.get('CATALOGUE', 'fsource_img', fallback=None), str)
    marzred          = _parse_val(config.get('CATALOGUE', 'marzred', fallback=None), str)
    delta_eo_sncut   = _parse_val(config.get('CATALOGUE', 'delta_eo_sncut', fallback='0.5,0.5'), 'tuple_float')
    sn_eo_cut        = _parse_val(config.get('CATALOGUE', 'sn_eo_cut', fallback='3.0,3.0'), 'tuple_float')
    fracnpix         = _parse_val(config.get('CATALOGUE', 'fracnpix', fallback=None), float)

    # =========================================================================
    # STEP 1 — Extraction
    # =========================================================================
    products = extract(
        fcube_input=fcube,
        fvar_input=fvar,
        extdata=extdata,
        extvar=extvar,
        mask2d=mask2d,
        mask2dpost=mask2dpost,
        snthreshold=snthreshold,
        spatsmooth=spatsmooth,
        specsmooth=specsmooth,
        connectivity=connectivity,
        maskspedge=maskspedge,
        mindz=mindz,
        maxdz=maxdz,
        minvox=minvox,
        minarea=minarea,
        zmin=zmin,
        zmax=zmax,
        lmin=lmin,
        lmax=lmax,
        outdir=outdir_extraction,
        do_continuum_sub=do_continuum_sub,
        rebinfac=rebinfac,
        filtsize=filtsize
    )

    # =========================================================================
    # STEP 2 — Covariance Estimation
    # =========================================================================
    # Use the labels map produced in Step 1 as a source mask if no custom one was provided
    actual_mask_cov = mask_source_cov if mask_source_cov is not None else products['fsegmap']

    print('\n' + '=' * 60)
    print('Empirical noise-covariance estimation')
    print('=' * 60)

    covmap, ssegments, wsegments, fitres = covariance(
        fcube=products['fcube_filtered'],
        fvar=products['fvar_filtered'],
        outdir=outdir_covariance,
        extcube=0,
        extvar=0,
        allsizes=allsizes,
        nsegments=nsegments,
        dl=dl,
        nsamples=nsamples,
        mask_source=actual_mask_cov,
        pixel_scale=pixel_scale,
        fitdeg=fitdeg,
        plot=plot_covariance,
        verbose=verbose_cov
    )

    # =========================================================================
    # STEP 3 — Catalogue Building
    # =========================================================================
    print('\n' + '=' * 60)
    print('Building the final catalogue')
    print('=' * 60)

    final_cat = build_em_catalog(
        fcube_for_extraction=products['fcube_filtered'],
        fvar_for_extraction=products['fvar_filtered'],
        fsegmap=products['fsegmap'],
        catpath=products['fcatalogue'],
        outdir=outdir_catalogue,
        cov_dir=outdir_covariance,
        target_z=target_z,
        rest_line=rest_line,
        vel_cut=vel_cut,
        fcube_median=fcube_median,
        fcube_odd=fcube_odd,
        fcube_even=fcube_even,
        fcube_median_var=fcube_median_var,
        fcube_odd_var=fcube_odd_var,
        fcube_even_var=fcube_even_var,
        fcube_for_spectra=fcube_for_spectra,
        fsource_img=fsource_img,
        marzred=marzred,
        SNcut=sncut,
        DeltaEOSNcut=delta_eo_sncut,
        SNEOcut=sn_eo_cut,
        fracnpix=fracnpix,
        checkimg=checkimg,
        padding=padding,
        pixel_scale=pixel_scale
    )

    print('\n' + '=' * 60)
    print('PIPELINE DONE')
    print('=' * 60)
    print(f"Total emitters in selected catalogue: {len(final_cat)}")
    print(f"Output saved to: {outdir_catalogue}")


if __name__ == '__main__':
    main()
