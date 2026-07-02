#!/usr/bin/env python
# coding: utf-8
# AUTHORS: MF, DT
# VERSION: 2.0
#
# Find_Em_SHINE/catalogue.py
# ---------------------------
# Final emitter catalogue construction, S/N correction, source images,
# and spectral extraction.

import os
import warnings
from pathlib import Path

import numpy as np

from astropy.io import fits, ascii
from astropy.table import Column, Table

# mypython is required for spectral extraction inside build_emitter_catalogue.
# The import is deferred to runtime so that the package can be loaded
# without mypython installed (spectral extraction is simply skipped).
_utl = None
try:
    from mypython.ifu import muse_utils as _utl
except ImportError:
    pass


warnings.filterwarnings("ignore", category=RuntimeWarning)


# =============================================================================
# 1.  velocityoffset — velocity offset from a target redshift
# =============================================================================

def velocityoffset(lambda_obs, target_z, rest_line):
    """Compute the line-of-sight velocity offset relative to a target redshift.

    For each observed wavelength ``lambda_obs``, the function assumes the
    emission line is ``rest_line`` and computes the redshift of the source.
    The velocity offset is then measured relative to ``target_z``.

    Parameters
    ----------
    lambda_obs : array-like
        Observed wavelength(s) in Å (e.g., the ``LambdaFL`` column of a SHINE
        catalogue).
    target_z : float
        Reference (target) redshift.
    rest_line : float
        Rest-frame wavelength of the emission line in Å.

    Returns
    -------
    veloffset : numpy.ndarray
        Line-of-sight velocity offset in km/s.  Positive values correspond
        to sources receding faster than the target (redshifted).

    Notes
    -----
    The formula used is the non-relativistic Doppler approximation:

    .. math::

        z_\\text{src} = \\frac{\\lambda_\\text{obs}}{\\lambda_\\text{rest}} - 1

        v = c \\cdot \\frac{z_\\text{src} - z_\\text{target}}{1 + z_\\text{target}}
    """
    c_kms    = 299792.458  # km/s
    z_src    = np.asarray(lambda_obs) / rest_line - 1.0
    veloffset = c_kms * (z_src - target_z) / (1.0 + target_z)
    return veloffset


# =============================================================================
# 2.  compute_corrected_snr — per-source S/N corrected for correlated noise
# =============================================================================

def compute_corrected_snr(catalog, covariance, segmap, cube, var,
                        cube_med=None, cube_odd=None, cube_even=None, 
                        var_med=None, var_odd=None, var_even=None,
                        apermap=None, contzcat=None):
                        
    """Compute the effective S/N of SHINE detections corrected for correlated noise.

    For each source in the catalogue the raw flux S/N is measured directly
    on the voxels belonging to that source in the segmentation map, and then
    divided by the covariance correction factor (which accounts for the spatial
    correlation of the noise introduced by the seeing or spatial smoothing).

    Optional half-exposure cubes (``cube_odd``, ``cube_even``) and a median
    cube (``cube_med``) can be provided so that the same S/N metric is computed
    on those products.  An aperture map with associated redshift catalogue can
    be used to flag sources that overlap with known continuum emitters.

    Parameters
    ----------
    catalog : astropy.table.Table
        Source catalogue produced by SHINE, with columns ``ID``, ``Xcent``
        (or ``XcentFL``), ``Ycent`` (or ``YcentFL``), ``Zcent``
        (or ``ZcentFL``), ``Xmin``, ``Xmax``, ``Ymin``, ``Ymax``,
        ``Zmin``, ``Zmax``.
        The catalogue must already contain the columns ``SNR``, ``SNR_odd``,
        ``SNR_even``, ``SNR_med``, ``covfac``, and ``BoxFraction``
        (initialised to zero) and ``OverContinuum`` (initialised to False).
    covariance : array-like, shape (N,)
        Per-source covariance correction factors. An array of ones disables
        the correction (pure Poisson noise).
    segmap : numpy.ndarray, shape (nz, ny, nx)
        3-D segmentation cube produced by SHINE, where each voxel contains
        the ID of the source it belongs to (0 = background).
    cube : numpy.ndarray, shape (nz, ny, nx)
        Main data cube used for extraction.
    var : numpy.ndarray, shape (nz, ny, nx)
        Variance of the main data cube.
    cube_med : numpy.ndarray or None, optional
        Median data cube. If provided, ``SNR_med`` is filled. Default *None*.
    cube_odd : numpy.ndarray or None, optional
        Odd-exposure data cube. If provided, ``SNR_odd`` is filled.
        Default *None*.
    cube_even : numpy.ndarray or None, optional
        Even-exposure data cube. If provided, ``SNR_even`` is filled.
        Default *None*.
    var_med : numpy.ndarray or None, optional
        Variance for the median cube. Required if ``cube_med`` is provided.
    var_odd : numpy.ndarray or None, optional
        Variance for the odd cube. Required if ``cube_odd`` is provided.
    var_even : numpy.ndarray or None, optional
        Variance for the even cube. Required if ``cube_even`` is provided.
    apermap : numpy.ndarray or None, optional
        2-D aperture map of continuum-detected sources (integer IDs). If
        provided together with ``contzcat``, the ``OverContinuum`` flag is set
        for sources whose centroid falls on a continuum source with a reliable
        redshift (quality flag > 2 in ``contzcat``).
    contzcat : astropy.table.Table or None, optional
        Redshift catalogue for continuum sources in Marz format, with columns
        ``#ID`` and ``QOP``. Required if ``apermap`` is provided.

    Returns
    -------
    catalog : astropy.table.Table
        The input catalogue with updated columns ``SNR``, ``covfac``,
        ``BoxFraction``, ``OverContinuum``, and (if the corresponding cubes
        are provided) ``SNR_med``, ``SNR_odd``, ``SNR_even``.
    """
    for j in range(len(catalog)):

        Id = catalog['ID'][j]

        # Flux-weighted centroid (fall back to geometric if not available)
        try:
            x = int(catalog['XcentFL'][j])
        except KeyError:
            x = int(catalog['Xcent'][j])
        try:
            y = int(catalog['YcentFL'][j])
        except KeyError:
            y = int(catalog['Ycent'][j])
        try:
            z = int(catalog['ZcentFL'][j])
        except KeyError:
            z = int(catalog['Zcent'][j])

        x1 = int(catalog['Xmin'][j])
        x2 = int(catalog['Xmax'][j])
        y1 = int(catalog['Ymin'][j])
        y2 = int(catalog['Ymax'][j])
        z1 = int(catalog['Zmin'][j])
        z2 = int(catalog['Zmax'][j])

        cutsegmap = segmap[z1:z2, y1:y2, x1:x2]
        okpix     = (cutsegmap == Id)

        # Fraction of voxels within a 5×5×5 bounding box
        cutcutsegmap = cutsegmap[
            z - 2 - z1: z + 3 - z1,
            y - 2 - y1: y + 3 - y1,
            x - 2 - x1: x + 3 - x1,
        ]
        cutokpix = (cutcutsegmap == Id)
        total_okpix = np.sum(okpix)
        catalog['BoxFraction'][j] = (
            float(np.sum(cutokpix)) / total_okpix
            if total_okpix > 0 else 0.0
        )

        # Main S/N corrected for covariance
        raw_snr           = (
            np.nansum(cube[z1:z2, y1:y2, x1:x2][okpix])
            / np.sqrt(np.nansum(var[z1:z2, y1:y2, x1:x2][okpix]))
        )
        catalog['SNR'][j]    = raw_snr / covariance[j]
        catalog['covfac'][j] = covariance[j]

        # Median cube S/N
        if cube_med is not None:
            if var_med is None:
                raise ValueError(
                    'var_med must be provided when cube_med is set.'
                )
            snr_med = (
                np.nansum(cube_med[z1:z2, y1:y2, x1:x2][okpix])
                / np.sqrt(np.nansum(var_med[z1:z2, y1:y2, x1:x2][okpix]))
            )
            catalog['SNR_med'][j] = snr_med / covariance[j]

        # Odd cube S/N
        if cube_odd is not None:
            if var_odd is None:
                raise ValueError(
                    'var_odd must be provided when cube_odd is set.'
                )
            snr_odd = (
                np.nansum(cube_odd[z1:z2, y1:y2, x1:x2][okpix])
                / np.sqrt(np.nansum(var_odd[z1:z2, y1:y2, x1:x2][okpix]))
            )
            catalog['SNR_odd'][j] = snr_odd / covariance[j]

        # Even cube S/N
        if cube_even is not None:
            if var_even is None:
                raise ValueError(
                    'var_even must be provided when cube_even is set.'
                )
            snr_even = (
                np.nansum(cube_even[z1:z2, y1:y2, x1:x2][okpix])
                / np.sqrt(np.nansum(var_even[z1:z2, y1:y2, x1:x2][okpix]))
            )
            catalog['SNR_even'][j] = snr_even / covariance[j]

        # Continuum overlap flag
        if apermap is not None and contzcat is not None:
            cntid = apermap[y, x]
            if cntid > 0:
                zid = (contzcat['#ID'] == int(cntid))
                if contzcat['QOP'][zid] > 2:
                    catalog['OverContinuum'][j] = True

    return catalog


# =============================================================================
# 3.  make_source_images — per-source cutout images from the segmentation map
# =============================================================================

def make_source_images(cubelist, segcube, header, catentry, Id, outdir,
                       outnamelist, padding=0):
    """Generate FITS cutout images for a single detected source.

    Iterates over a list of data cubes (e.g., mean, median, odd, even), and
    for each cube produces a 2-D image by collapsing the voxels associated
    with the given source ID in the segmentation map.  The result is stored
    as a multi-extension FITS file.

    Parameters
    ----------
    cubelist : list of numpy.ndarray
        List of 3-D data cubes (shape ``(nz, ny, nx)``).  All cubes must share
        the same spatial and spectral dimensions as ``segcube``.
    segcube : numpy.ndarray, shape (nz, ny, nx)
        Segmentation cube produced by SHINE.  Each voxel contains the integer
        ID of the associated source (0 = background).
    header : astropy.io.fits.Header
        FITS header of the data cube (used to propagate WCS keywords).
    catentry : astropy.table.Row
        A single row from the SHINE output catalogue, containing at least
        the columns ``ZcentFL`` (or ``Zcent``), ``XcentFL`` (or ``Xcent``),
        ``YcentFL`` (or ``Ycent``), ``Xmin``, ``Xmax``, ``Ymin``, ``Ymax``,
        ``Zmin``, ``Zmax``.
    Id : int
        Source identifier to process.
    outdir : str
        Directory where the output FITS file is written.
    outnamelist : list of str
        List of name tags, one per cube in ``cubelist``, used to label each
        image extension (e.g., ``['_mean', '_median', '_half1', '_half2']``).
    padding : int, optional
        Number of pixels to extend the spatial bounding box on each side to
        produce postage-stamp images.  Set to 0 (default) to skip postage
        stamps.

    Returns
    -------
    None
        Writes ``{outdir}/id{Id}_img.fits``.

    Notes
    -----
    The output FITS file has the following extension layout:

    * Extension 0 (primary): empty
    * Extension 1: ``segcube`` sub-cube trimmed to the padded bounding box
    * Extension 2: segmentation image (2-D collapse of the source mask)
    * Extensions 3…N+2: one 2-D image per cube in ``cubelist``
    * Extensions N+3…2N+3 (only if ``padding > 0``): padded postage stamps
    """
    mz, my, mx = segcube.shape
    pixmask = (segcube == Id) * 1.0

    # Flux-weighted centroid (fall back to geometric)
    try:
        zgeo = int(catentry['ZcentFL'])
    except KeyError:
        zgeo = int(catentry['Zcent'])
    try:
        xgeo = int(catentry['XcentFL'])
    except KeyError:
        xgeo = int(catentry['Xcent'])
    try:
        ygeo = int(catentry['YcentFL'])
    except KeyError:
        ygeo = int(catentry['Ycent'])

    x1 = int(catentry['Xmin'])
    x2 = int(catentry['Xmax'])
    y1 = int(catentry['Ymin'])
    y2 = int(catentry['Ymax'])
    z1 = int(catentry['Zmin'])
    z2 = int(catentry['Zmax'])

    # Padded bounding box
    xpad1 = max(x1 - padding, 0)
    xpad2 = min(x2 + padding, mx)
    ypad1 = max(y1 - padding, 0)
    ypad2 = min(y2 + padding, my)
    zpad1 = max(z1 - 2, 0)
    zpad2 = min(z2 + 2, mz)

    # Build sub-headers (remove spectral axis for 2-D products)
    imahdr = header.copy()
    for key in list(imahdr.keys()):
        if key.endswith('3') and key.startswith('C'):
            del imahdr[key]

    imapadhdr = imahdr.copy()
    imapadhdr['CRPIX1'] = imahdr.get('CRPIX1', 1) - xpad1
    imapadhdr['CRPIX2'] = imahdr.get('CRPIX2', 1) - ypad1

    cubpadhdr = header.copy()
    cubpadhdr['CRPIX1'] = header.get('CRPIX1', 1) - xpad1
    cubpadhdr['CRPIX2'] = header.get('CRPIX2', 1) - ypad1
    cubpadhdr['CRPIX3'] = header.get('CRPIX3', 1) - zpad1

    # Primary HDU (empty)
    hdu_zero = fits.PrimaryHDU([])

    # Segmentation sub-cube
    segmapshort = pixmask[zpad1:zpad2, ypad1:ypad2, xpad1:xpad2]
    hdu_segcube = fits.ImageHDU(segmapshort, header=cubpadhdr)

    # Segmentation image (2-D projection)
    segima      = np.nansum(pixmask, axis=0)
    hdu_img_det = fits.ImageHDU(segima, header=imahdr)

    hdu_img_list = []
    hdu_pad_list = []

    for ii, thiscube in enumerate(cubelist):
        # Collapse source voxels; fill zeros with the centroid layer
        thisima        = np.nansum(thiscube * pixmask, axis=0)
        empty          = (thisima == 0)
        thisima[empty] = thiscube[zgeo, empty]

        hdu_img = fits.ImageHDU(thisima, header=imahdr)
        hdu_img_list.append(hdu_img)

        if padding > 0:
            hdu_pad = fits.ImageHDU(
                thisima[ypad1:ypad2, xpad1:xpad2], header=imapadhdr
            )
            hdu_pad_list.append(hdu_pad)

    # Assemble and write
    all_hdus = [hdu_zero, hdu_segcube, hdu_img_det] + hdu_img_list
    if padding > 0:
        pad_segima  = fits.ImageHDU(segima[ypad1:ypad2, xpad1:xpad2],
                                    header=imapadhdr)
        all_hdus   += [pad_segima] + hdu_pad_list

    fits.HDUList(all_hdus).writeto(
        os.path.join(outdir, f'id{Id}_img.fits'), overwrite=True
    )


# =============================================================================
# 4.  build_emitter_catalogue — full post-processing pipeline
# =============================================================================

def build_emitter_catalogue(
    fcube,
    fcube_var,
    fsegmap,
    catpath,
    outdir='./',
    cov_poly=None,
    cov_dir=None,
    target_z=None,
    rest_line=None,
    vel_cut=None,
    fcube_median=None,
    fcube_odd=None,
    fcube_even=None,
    fcube_median_var=None,
    fcube_odd_var=None,
    fcube_even_var=None,
    fcube_orig=None,
    fsource_img=None,
    marzred=None,
    SNcut=(7, 5),
    DeltaEOSNcut=(0.5, 0.5),
    SNEOcut=(3, 3),
    fracnpix=None,
    derived=True,
    checkimg=True,
    mask=None,
    startind=0,
    padding=50,
    pixel_scale=0.2,
):
    """Build the final catalogue of line emitters from a SHINE extraction run.

    Starting from the raw SHINE catalogue and segmentation map, this function:

    1. Initialises quality-assessment columns (S/N, even/odd S/N, covariance
       factor, velocity offset, continuum overlap flag, …).
    2. Optionally computes the velocity offset of each source with respect to
       a target redshift (requires ``target_z`` and ``rest_line``).
    3. Computes per-source S/N accounting for correlated noise via the
       covariance correction factor supplied in ``cov_poly``.
    4. Applies S/N cuts and optional additional quality criteria to assign a
       confidence class to each source.
    5. Generates per-source image cutouts (requires ``checkimg=True``).
    6. Optionally extracts 1-D spectra via the ``mypython.ifu.muse_utils``
       ``cube2spec`` routine if the ``mypython`` package is available and
       ``fcube_orig`` is provided.

    Parameters
    ----------
    fcube : str
        Path to the filtered data cube (``*FILTER_out.fits`` from SHINE).
    fcube_var : str
        Path to the filtered variance cube.
    fsegmap : str
        Path to the SHINE segmentation map (``*LABELS_out.fits``).
    catpath : str
        Path to the SHINE output catalogue (``*CATALOGUE_out.fits``).
    outdir : str, optional
        Directory where output files are written. Default is ``'./'``.
    cov_poly : numpy.ndarray or None, optional
        Covariance correction model.  Two formats are accepted:

        * 1-D array ``(N2, N1, N0)``: a single polynomial in aperture size
          (in arcseconds) applied to all sources.
        * 2-D array, shape ``(K, fitdeg + 3)``: one polynomial per wavelength
          interval.  Each row contains ``[Wmin, Wmax, P_n, …, P_0]`` as
          written by :func:`estimate_empirical_covariance`.

        If *None*, a correction of 1 (no correction) is applied.
    cov_dir : str or None, optional
        Path to the covariance estimation directory. If provided and ``cov_poly``
        is ``None``, the 1-D polynomial fit coefficients are automatically loaded
        from ``{cov_dir}/covariance_1Dfit.txt``. Default is ``None``.
    target_z : float or None, optional
        Target redshift for velocity-offset computation. Requires
        ``rest_line``.
    rest_line : float or None, optional
        Rest-frame wavelength of the emission line in Å.
    vel_cut : float or None, optional
        If set (together with ``target_z`` and ``rest_line``), removes sources
        with ``|veloffset| > vel_cut`` (km/s).
    fcube_median : str or None, optional
        Path to the median data cube. Enables ``SNR_med`` computation.
    fcube_odd : str or None, optional
        Path to the odd-exposure data cube. Enables ``SNR_odd`` computation.
    fcube_even : str or None, optional
        Path to the even-exposure data cube. Enables ``SNR_even`` computation.
    fcube_median_var : str or None, optional
        Variance for the median cube.
    fcube_odd_var : str or None, optional
        Variance for the odd cube.
    fcube_even_var : str or None, optional
        Variance for the even cube.
    fcube_orig : str or None, optional
        Path to the *unfiltered* data cube. If provided and ``mypython`` is
        installed, 1-D spectra are extracted for each source.
    fsource_img : str or None, optional
        Path to a 2-D aperture map of continuum sources (integer IDs). Used
        together with ``marzred`` to set the ``OverContinuum`` flag.
    marzred : str or None, optional
        Path to a Marz-format redshift catalogue for continuum sources.
    SNcut : tuple of float, optional
        S/N thresholds defining confidence classes. Sources are assigned to
        the first class whose threshold they meet. Default is ``(7, 5)``.
    DeltaEOSNcut : tuple of float, optional
        Maximum fractional even/odd S/N difference allowed per class. Only
        used when ``fcube_odd`` and ``fcube_even`` are provided. Default is
        ``(0.5, 0.5)``.
    SNEOcut : tuple of float, optional
        Minimum even and odd S/N per class. Default is ``(3, 3)``.
    fracnpix : float or None, optional
        If set, removes sources with fewer than ``fracnpix`` of their voxels
        within a 5×5×5 bounding box centred on the centroid.
    derived : bool, optional
        If *True* (default), compute S/N and other derived quantities.
    checkimg : bool, optional
        If *True* (default), generate per-source image cutouts.
    mask : str or None, optional
        Path to a 2-D FITS mask (1 = bad). Sources whose centroid falls on
        a masked pixel are removed from the catalogue.
    startind : int, optional
        Index from which to start the image-cutout loop (useful for
        restarting interrupted runs). Default is 0.
    padding : int, optional
        Spatial padding (pixels) for image cutouts. Default is 50.
    pixel_scale : float, optional
        Pixel scale in arcseconds per pixel, used to compute projected
        aperture sizes for the covariance correction. Default is 0.2.

    Returns
    -------
    catalog : astropy.table.Table
        Final catalogue with all quality-assessment columns and confidence
        class assignments.

    Notes
    -----
    Two catalogue files are written to ``outdir``:

    * ``{catname}_all_SNR.fits``    — all sources above the minimum S/N cut
      with derived metrics.
    * ``{catname}_select_SNR.fits`` — sources after all quality cuts, with
      confidence class assignments.

    Per-source image cutouts are written to ``{outdir}/objs/id{Id}/``.
    """
    catname = Path(catpath).stem
    cat_name = os.path.join(
        outdir, catname.split('.fits')[0] + '_select_SNR.fits'
    )

    os.makedirs(outdir, exist_ok=True)

    # ------------------------------------------------------------------
    # Load or recompute catalogue with derived quantities
    # ------------------------------------------------------------------
    if os.path.isfile(cat_name):
        print(f'Found existing selected catalogue: {cat_name}. '
              'Loading from disk (skipping derivation step).')
        catalog = Table.read(cat_name)
    else:
        # Load the raw SHINE catalogue
        print(f'Reading SHINE catalogue from: {catpath}')
        catalog = Table.read(catpath, format='fits')

        # Initialise quality columns
        n = len(catalog)
        catalog.add_columns([
            Column(np.zeros(n, dtype=float), name='SNR'),
            Column(np.zeros(n, dtype=float), name='SNR_odd'),
            Column(np.zeros(n, dtype=float), name='SNR_even'),
            Column(np.zeros(n, dtype=float), name='SNR_med'),
            Column(np.zeros(n, dtype=float), name='covfac'),
            Column(np.zeros(n, dtype=float), name='confidence'),
            Column(np.zeros(n, dtype=float), name='veloffset'),
            Column(np.zeros(n, dtype=float), name='EODeltaSN'),
            Column(np.zeros(n, dtype=float), name='BoxFraction'),
            Column(np.zeros(n, dtype=bool),  name='OverContinuum'),
        ])

        # ------------------------------------------------------------------
        # Covariance correction vector
        # ------------------------------------------------------------------
        if cov_poly is None and cov_dir is not None:
            fitcoeffs_path = os.path.join(cov_dir, 'covariance_1Dfit.txt')
            if os.path.isfile(fitcoeffs_path):
                print(f'Loading covariance polynomial fit from: {fitcoeffs_path}')
                cov_poly = np.loadtxt(fitcoeffs_path)
            else:
                warnings.warn(
                    f'cov_dir was provided, but covariance_1Dfit.txt '
                    f'was not found in {cov_dir}. Proceeding with no covariance correction.',
                    UserWarning
                )

        if cov_poly is None:
            covariance = np.ones(n, dtype=float)


        elif np.ndim(cov_poly) == 1:
            # Single polynomial in projected size (arcsec)
            size       = np.sqrt(catalog['Nspat']) * pixel_scale
            covariance = np.polyval(cov_poly, size)

        elif np.ndim(cov_poly) == 2:
            # Per-wavelength-interval polynomial
            # Each row: [Wmin, Wmax, P_n, ..., P_0]
            size       = np.sqrt(catalog['Nspat']) * pixel_scale
            covariance = np.ones(n, dtype=float)
            for ii in range(n):
                try:
                    lam = catalog['LambdaFL'][ii]
                except KeyError:
                    lam = catalog['Lambda'][ii]
                # Find the wavelength interval
                okind = np.where(lam > cov_poly[:, 0])[0]
                if len(okind) == 0:
                    okind = 0
                else:
                    okind = okind[-1]
                # Polynomial coefficients start at column 2
                covariance[ii] = np.polyval(cov_poly[okind, 2:], size[ii])
        else:
            raise ValueError(
                'cov_poly must be None, a 1-D array, or a 2-D array.'
            )

        # ------------------------------------------------------------------
        # Open cubes
        # ------------------------------------------------------------------
        def _open_cube(fpath, label=''):
            if fpath is None:
                return None
            hdu = fits.open(fpath)
            try:
                data = hdu[1].data
            except IndexError:
                data = hdu[0].data
            hdu.close()
            print(f'  -> {label} in {fpath}')
            return data

        print('Reading cubes …')
        cubehdu = fits.open(fcube)
        try:
            cube    = cubehdu[1].data
            cubehdr = cubehdu[1].header
        except IndexError:
            cube    = cubehdu[0].data
            cubehdr = cubehdu[0].header

        cube_var    = _open_cube(fcube_var,    'filtered variance')
        segmap      = fits.open(fsegmap)[0].data
        cube_odd    = _open_cube(fcube_odd,    'odd cube')
        cube_odd_v  = _open_cube(fcube_odd_var, 'odd variance')
        cube_even   = _open_cube(fcube_even,   'even cube')
        cube_even_v = _open_cube(fcube_even_var, 'even variance')
        cube_median = _open_cube(fcube_median, 'median cube')
        cube_med_v  = _open_cube(fcube_median_var, 'median variance')

        apermap  = None
        contzcat = None
        if fsource_img is not None and marzred is not None:
            apermap = fits.open(fsource_img)[0].data
            try:
                contzcat = ascii.read(marzred, format='csv', header_start=2)
            except Exception as exc:
                raise ValueError(
                    f'Could not read the Marz redshift file {marzred}: {exc}'
                )

        # ------------------------------------------------------------------
        # Velocity offset + optional trim
        # ------------------------------------------------------------------
        if target_z is not None and rest_line is not None:
            try:
                lam_col = catalog['LambdaFL']
            except KeyError:
                lam_col = catalog['Lambda']
            catalog['veloffset'] = velocityoffset(lam_col, target_z, rest_line)

            if vel_cut is not None:
                select  = np.abs(catalog['veloffset']) <= vel_cut
                catalog = catalog[select]
                # Re-slice covariance to match trimmed catalogue
                covariance = covariance[select]

        # ------------------------------------------------------------------
        # Derived S/N and metrics
        # ------------------------------------------------------------------
        if derived:
            print(f'\nCalculating corrected S/N for {len(catalog)} sources …')
            catalog = compute_corrected_snr(
                catalog, covariance, segmap, cube, cube_var,
                cube_med=cube_median, cube_odd=cube_odd, cube_even=cube_even,
                var_med=cube_med_v, var_odd=cube_odd_v, var_even=cube_even_v,
                apermap=apermap, contzcat=contzcat,
            )

        # Even/odd S/N fractional difference
        if cube_even is not None and cube_odd is not None:
            rel_diff = (
                np.abs(catalog['SNR_even'] - catalog['SNR_odd'])
                / np.minimum(catalog['SNR_even'], catalog['SNR_odd'])
            )
            catalog['EODeltaSN'] = rel_diff

        # Write full catalogue (before cuts)
        all_cat_path = os.path.join(
            outdir, catname.split('.fits')[0] + '_all_SNR.fits'
        )
        print(f'\nWriting full catalogue to {all_cat_path}')
        catalog.write(all_cat_path, format='fits', overwrite=True)

        # ------------------------------------------------------------------
        # Minimum S/N cut
        # ------------------------------------------------------------------
        select  = catalog['SNR'] >= np.amin(SNcut)
        catalog = catalog[select]
        print(f'Sources above minimum S/N ({np.amin(SNcut)}): {len(catalog)}')

        # Optional spatial mask
        if mask is not None:
            hdu_mask = fits.open(mask)
            try:
                msk = hdu_mask[0].data
            except IndexError:
                msk = hdu_mask[1].data
            hdu_mask.close()
            try:
                masked  = msk[
                    np.array(catalog['YcentFL'], dtype=int),
                    np.array(catalog['XcentFL'], dtype=int),
                ]
            except KeyError:
                masked  = msk[
                    np.array(catalog['Ycent'], dtype=int),
                    np.array(catalog['Xcent'], dtype=int),
                ]
            catalog = catalog[masked == 0]

        # ------------------------------------------------------------------
        # Confidence class assignment
        # ------------------------------------------------------------------
        print('\nAssigning confidence classes …')
        for iclass, iSN in enumerate(SNcut):

            have_eo = cube_even is not None and cube_odd is not None
            have_fp = fracnpix is not None
            have_ct = fsource_img is not None

            if have_eo and have_fp and have_ct:
                sel = (
                    (catalog['SNR']          >= iSN)
                    & (catalog['SNR_odd']    >= SNEOcut[iclass])
                    & (catalog['SNR_even']   >= SNEOcut[iclass])
                    & (catalog['EODeltaSN']  <= DeltaEOSNcut[iclass])
                    & (~catalog['OverContinuum'])
                    & (catalog['BoxFraction'] >= fracnpix)
                    & (catalog['confidence'] == 0)
                )
            elif have_eo and have_fp:
                sel = (
                    (catalog['SNR']          >= iSN)
                    & (catalog['SNR_odd']    >= SNEOcut[iclass])
                    & (catalog['SNR_even']   >= SNEOcut[iclass])
                    & (catalog['EODeltaSN']  <= DeltaEOSNcut[iclass])
                    & (catalog['BoxFraction'] >= fracnpix)
                    & (catalog['confidence'] == 0)
                )
            elif have_eo:
                sel = (
                    (catalog['SNR']         >= iSN)
                    & (catalog['SNR_odd']   >= SNEOcut[iclass])
                    & (catalog['SNR_even']  >= SNEOcut[iclass])
                    & (catalog['EODeltaSN'] <= DeltaEOSNcut[iclass])
                    & (catalog['confidence'] == 0)
                )
            else:
                sel = (
                    (catalog['SNR'] >= iSN)
                    & (catalog['confidence'] == 0)
                )

            catalog['confidence'][sel] = iclass + 1

        # Write selected catalogue
        print(f'Writing selected catalogue to {cat_name}')
        catalog.write(cat_name, format='fits', overwrite=True)

        cubehdu.close()

    # ------------------------------------------------------------------
    # Per-source image cutouts
    # ------------------------------------------------------------------
    if checkimg:
        objs_dir = os.path.join(outdir, 'objs')
        os.makedirs(objs_dir, exist_ok=True)

        # Re-open cubes needed for image generation
        cubehdu = fits.open(fcube)
        try:
            cube    = cubehdu[1].data
            cubehdr = cubehdu[1].header
        except IndexError:
            cube    = cubehdu[0].data
            cubehdr = cubehdu[0].header

        segmap = fits.open(fsegmap)[0].data

        cube_median = (
            fits.open(fcube_median)[1 if fits.open(fcube_median)[1:] else 0].data
            if fcube_median else None
        )
        cube_odd = (
            fits.open(fcube_odd)[1 if fits.open(fcube_odd)[1:] else 0].data
            if fcube_odd else None
        )
        cube_even = (
            fits.open(fcube_even)[1 if fits.open(fcube_even)[1:] else 0].data
            if fcube_even else None
        )

        if cube_median is not None and cube_odd is not None and cube_even is not None:
            hdulist     = [cube, cube_median, cube_odd, cube_even]
            outnamelist = ['_mean', '_median', '_half1', '_half2']
        else:
            hdulist     = [cube]
            outnamelist = ['_mean']

        total = len(catalog)
        step  = max(total // 10, 1)

        if fcube_orig is not None:
            print(f'\nExtracting images and spectra for {total} sources')
        else:
            print(f'\nExtracting images for {total} sources')

        for ii in range(startind, total):
            objid  = catalog['ID'][ii]
            objdir = os.path.join(objs_dir, f'id{objid}')

            if os.path.isdir(objdir):
                print(f'Output for id{objid} already exists — skipping.')
                continue

            os.makedirs(objdir, exist_ok=True)

            make_source_images(
                hdulist, segmap, cubehdr, catalog[ii],
                objid, objdir + os.sep, outnamelist, padding=padding,
            )

            # Spectral extraction via mypython.cube2spec
            if fcube_orig is not None:
                if _utl is None:
                    raise ImportError(
                        'mypython is required for spectral extraction. '
                        'Install it with "pip install mypython" or set '
                        'fcube_orig=None to skip spectral extraction.'
                    )
                savename = os.path.join(objdir, 'spectrum.fits')
                _utl.cube2spec(
                    fcube_orig, 0.0, 0.0, 0.0,
                    shape='mask', helio=0, mask=segmap,
                    twod=True, tovac=True, write=savename,
                    idsource=objid,
                )

            if (ii - startind) % step == 0:
                progress = (ii - startind) / total * 100
                print(f'Progress: {progress:.0f}%  (source {ii + 1}/{total})')

        cubehdu.close()

    return catalog
