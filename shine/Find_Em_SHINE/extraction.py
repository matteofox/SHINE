#!/usr/bin/env python
# coding: utf-8
# AUTHORS: MF, DT
# VERSION: 2.0
#
# Find_Em_SHINE/extraction.py
# ----------------------------
# Convenience wrapper around SHINE.runextraction() with defaults
# optimised for line-emitter detection.

import os
from pathlib import Path

import numpy as np
from astropy.io import fits

from ..SHINE import runextraction
from ..shine_utils import clean_clube


def extract(fcube, fvar, extdata=0, extvar=0,
            mask2d=None, mask2dpost=None,
            snthreshold=2.0, spatsmooth=2.0, specsmooth=0.0,
            connectivity=26, maskspedge=20,
            mindz=1, maxdz=200, minvox=1, minarea=1,
            zmin=None, zmax=None, lmin=None, lmax=None,
            outdir='./',
            do_continuum_sub=False, rebinfac=40, filtsize=7):

    """Run SHINE extraction with defaults optimised for line-emitter detection.

    This is a convenience wrapper around :func:`shine.SHINE.runextraction`
    that:

    * optionally performs continuum subtraction before extraction
      (via :func:`shine.shine_utils.clean_clube`);
    * forces output of labels, filtered data, and filtered variance cubes
      (required by the downstream covariance and catalogue steps).

    Parameters
    ----------
    fcube : str
        Path to the science data cube FITS file.
    fvar : str
        Path to the variance cube FITS file, or a numeric string /
        ``'-1'`` (see :func:`shine.SHINE.runextraction`).
    extdata : int, optional
        HDU extension index for the data cube. Default is 0.
    extvar : int, optional
        HDU extension index for the variance cube. Default is 0.
    mask2d : str or None, optional
        Path to a 2-D pre-smoothing mask (1 = bad).
    mask2dpost : str or None, optional
        Path to a 2-D post-smoothing mask (1 = bad).
    snthreshold : float, optional
        S/N threshold for voxel inclusion. Default is 2.0.
    spatsmooth : float, optional
        Spatial Gaussian sigma (pixels). Default is 2.0.
    specsmooth : float, optional
        Spectral Gaussian sigma (pixels). Default is 0.0 (disabled).
    connectivity : int, optional
        Voxel connectivity: 6, 18, or 26. Default is 26.
    maskspedge : int, optional
        Pixels to mask around the field edges. Default is 20.
    mindz : int, optional
        Minimum spectral extent per source. Default is 1.
    maxdz : int, optional
        Maximum spectral extent per source. Default is 200.
    minvox : int, optional
        Minimum total voxels per source. Default is 1.
    minarea : int, optional
        Minimum projected spatial area per source. Default is 1.
    zmin : int or None, optional
        Starting layer index (0-based).
    zmax : int or None, optional
        Ending layer index (0-based).
    lmin : float or None, optional
        Starting wavelength (Å).
    lmax : float or None, optional
        Ending wavelength (Å).
    outdir : str, optional
        Output directory. Default is ``'./'``.
    do_continuum_sub : bool, optional
        If *True*, subtract the continuum before extraction using
        :func:`clean_clube`. Default is *False*.
    rebinfac : int, optional
        Spectral re-binning factor for continuum subtraction. Default is 40.
    filtsize : int, optional
        Median filter width for continuum subtraction. Default is 7.

    Returns
    -------
    products : dict
        Dictionary with paths to the SHINE output products:

        * ``'fcube_filtered'`` — filtered (smoothed) data cube
        * ``'fvar_filtered'``  — filtered (smoothed) variance cube
        * ``'fsegmap'``        — segmentation (labels) cube
        * ``'fcatalogue'``     — raw extraction catalogue
        * ``'fcube_clean'``    — continuum-subtracted cube
          (only if ``do_continuum_sub=True``, otherwise *None*)

    Examples
    --------
    >>> from shine.Find_Em_SHINE import extract
    >>> products = extract(
    ...     'Datacube.fits', 'Varcube.fits',
    ...     snthreshold=2.0, spatsmooth=2.0, outdir='./extraction/',
    ... )
    >>> print(products['fsegmap'])
    ./extraction/Datacube.LABELS_out.fits
    """
    os.makedirs(outdir, exist_ok=True)

    cube_stem = Path(fcube).stem
    var_stem  = Path(fvar).stem

    # ------------------------------------------------------------------
    # Optional continuum subtraction
    # ------------------------------------------------------------------
    fcube_clean = None
    if do_continuum_sub:
        print('\n' + '=' * 60)
        print('Continuum subtraction (clean_clube)')
        print('=' * 60)

        hducube = fits.open(fcube)
        cube    = hducube[extdata].data
        header  = hducube[extdata].header
        hducube.close()

        cube_clean = clean_clube(cube, filtsize=filtsize, rebinfac=rebinfac)

        fcube_clean = os.path.join(outdir, cube_stem + '_CLEAN.fits')
        hdu_out = fits.PrimaryHDU(np.array(cube_clean), header=header)
        hdu_out.header['HISTORY'] = (
            'Continuum subtracted with shine.shine_utils.clean_clube'
        )
        hdu_out.writeto(fcube_clean, overwrite=True)
        print(f'  Saved continuum-subtracted cube to: {fcube_clean}')

        fcube_for_extraction = fcube_clean
    else:
        fcube_for_extraction = fcube

    # ------------------------------------------------------------------
    # SHINE extraction
    # ------------------------------------------------------------------
    print('\n' + '=' * 60)
    print('SHINE extraction')
    print('=' * 60)
    print(f'  S/N threshold  : {snthreshold}')
    print(f'  Spatial smooth : σ = {spatsmooth} pix')
    print(f'  Spectral smooth: σ = {specsmooth} pix')
    print(f'  Connectivity   : {connectivity}')
    print(f'  Min voxels     : {minvox}')
    print(f'  Min Δz         : {mindz}')
    print(f'  Max Δz         : {maxdz}')
    print(f'  Min area       : {minarea}')
    print(f'  Edge mask      : {maskspedge} pix')
    print(f'  Output dir     : {outdir}')

    runextraction(
        fcube_for_extraction,
        fvar,
        extdata=extdata,
        extvardata=extvar,
        mask2d=mask2d,
        mask2dpost=mask2dpost,
        snthreshold=snthreshold,
        spatsmooth=spatsmooth,
        specsig=specsmooth,
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
        outdir=outdir,
        # Always write the outputs required by downstream steps
        writelabels=True,
        writesmdata=True,
        writesmvar=True,
        writesmsnr=False,
        writesubcube=False,
        writevardata=False,
    )

    # ------------------------------------------------------------------
    # Build output product paths
    # ------------------------------------------------------------------
    if do_continuum_sub:
        extr_stem = Path(fcube_for_extraction).stem
    else:
        extr_stem = cube_stem

    products = {
        'fcube_filtered': os.path.join(
            outdir, f'{extr_stem}.FILTER_out.fits'
        ),
        'fvar_filtered': os.path.join(
            outdir, f'{var_stem}.FILTER_out.fits'
        ),
        'fsegmap': os.path.join(
            outdir, f'{extr_stem}.LABELS_out.fits'
        ),
        'fcatalogue': os.path.join(
            outdir, f'{extr_stem}.CATALOGUE_out.fits'
        ),
        'fcube_clean': fcube_clean,
    }

    print(f'\n  Products written:')
    for key, val in products.items():
        if val is not None:
            print(f'    {key:20s}: {val}')

    return products
