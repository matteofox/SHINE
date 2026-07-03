#!/usr/bin/env python
# coding: utf-8
# AUTHORS: MF, DT
# VERSION: 2.0
#
# Find_Em_SHINE/covariance.py
# ----------------------------
# Empirical noise-covariance estimation for spectroscopic cubes.

import os
import warnings

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

# Optional dependency: sep is used only when no mask_source is provided
# (automatic source detection).
try:
    import sep
    _SEP_AVAILABLE = True
except ImportError:
    _SEP_AVAILABLE = False

warnings.filterwarnings("ignore", category=RuntimeWarning)


# =============================================================================
# estimate_empirical_covariance
# =============================================================================

def estimate_empirical_covariance(fcube_for_extraction, fvar_for_extraction, outdir, extcube=0, extvar=0, allsizes=None, nsegments=500, dl=4, nsamples=10000, max_attempts=100000, 
    mask_source=None, mask_edge=None, pixel_scale=0.2, fitdeg=2, plot=True, verbose=True, sampling=True, random_seed=0):
    
    """Estimate the empirical noise covariance of a MUSE-like spectroscopic cube.

    The function samples random apertures of increasing projected size at
    random positions and wavelengths throughout the cube. For each aperture
    it computes the measured flux divided by the expected standard deviation.
    (:math:`\\text{SNR} = F / \\sqrt{V}`). In the absence of correlated noise
    this ratio would be Gaussian with unit standard deviation; in practice it
    exceeds unity and depends on both the aperture size and the wavelength.
    This function measures that excess as a function of both dimensions and
    fits a per-wavelength polynomial model that can later be used to correct
    the S/N of extracted sources.

    Parameters
    ----------
    fcube_for_extraction : str
        Path to the FITS file containing the filtered data cube (ideally the
        ``*FILTER_out.fits`` product written by ``SHINE --writesmdata``).
    fvar_for_extraction : str
        Path to the FITS file containing the filtered variance cube (ideally
        the ``*FILTER_out.fits`` product written by ``SHINE --writesmvar``).
    outdir : str
        Directory where all output files will be written.  Created
        automatically if it does not exist.
    extcube : int, optional
        HDU extension index for the data cube. Default is 0.
    extvar : int, optional
        HDU extension index for the variance cube. Default is 0.
    allsizes : list of int, optional
        Aperture half-sizes (in pixels) to probe.  Each value defines a
        square aperture of ``size x size`` pixels.  If *None*, defaults to
        ``range(2, 31)`` (i.e., 2 to 30 pixels inclusive).
    nsegments : int, optional
        Number of wavelength segments into which the spectral range is divided.
        Default is 500.
    dl : int, optional
        Spectral depth (number of layers) of each sampled aperture. Default
        is 4.
    nsamples : int, optional
        Number of valid (unmasked) samples to collect per aperture size per
        wavelength segment before moving on.  Default is 10 000.
    max_attempts : int, optional
        Maximum number of random draws per (size, segment) combination before
        giving up. Default is 100 000.
    mask_source : str or None, optional
        Path to a FITS file containing a 2-D source-detection map (e.g., the
        ``*LABELS_out.fits`` output of a SHINE 2-D run). Non-zero pixels are
        masked.  If *None* and ``sep`` is installed, sources are detected
        automatically using :mod:`sep`. If *None* and ``sep`` is not installed,
        only the edge mask is applied.
    mask_edge : str or None, optional
        Path to a FITS file containing a 2-D edge mask (1 = bad, 0 = good).
        If *None*, the edge mask is derived automatically from NaN/zero pixels
        in the median-collapsed white-light image.
    pixel_scale : float, optional
        Pixel scale in arcseconds per pixel. Used to convert aperture sizes
        from pixels to arcseconds for the covariance fit and plots. Default
        is 0.2 (MUSE nominal).
    fitdeg : int, optional
        Degree of the polynomial fit to the covariance as a function of
        aperture size (in arcseconds) per wavelength segment. Default is 2.
    plot : bool, optional
        If *True*, save a 2-D colour map of the covariance as a function of
        aperture and wavelength. Default is *True*.
    verbose : bool, optional
        If *True*, print progress information. Default is *True*.
    sampling : bool, optional
        If *True* (default), run the full Monte Carlo sampling loop and
        save the results to the ``.npz`` file.  If *False*, skip the
        sampling entirely and reload ``savenormstd``, ``savestd``,
        ``allsizes``, and ``wsegments`` from the ``.npz`` file written by
        a previous run (the file must already exist in ``outdir``).
        This is useful when you only want to change the polynomial fit
        degree (``fitdeg``) or the plot without repeating the expensive
        sampling phase.
    random_seed : int or None, optional
        If set, fix the NumPy random seed before the sampling loop so that
        results are exactly reproducible across runs.  Default is 0.

    Returns
    -------
    covmap : numpy.ndarray, shape (nsegments, len(allsizes))
        Effective S/N correction factor (:math:`\\sigma_\\text{eff}/\\sigma_N`)
        per wavelength segment and aperture size.
    ssegments : numpy.ndarray
        Aperture sizes in arcseconds.
    wsegments : numpy.ndarray, shape (nsegments + 1,)
        Wavelength bin edges (Å).
    fitres : numpy.ndarray, shape (nsegments, fitdeg + 1)
        Polynomial coefficients of the 1-D fit to the covariance per segment.
        Each row corresponds to one wavelength segment; columns are the
        coefficients ``[P_n, ..., P_1, P_0]`` (highest degree first, as
        returned by :func:`numpy.polyfit`).

    Notes
    -----
    **Output files** written to ``outdir``:

    * ``npixscaling_spscale_{nsegments}seg.npz`` — compressed NumPy archive
      with keys ``nstd`` (covmap), ``std`` (raw flux dispersion per segment),
      ``size`` (allsizes array), ``waves`` (wsegments array).
    * ``covariance_map.pdf`` — 2-D colour map of the covariance
      (only if ``plot=True``).
    * ``covariance_1Dfit.txt`` — ASCII table with columns
      ``Wmin Wmax P2 P1 P0`` (the number of polynomial columns scales with
      ``fitdeg``).

    Examples
    --------
    >>> from shine.Find_Em_SHINE import covariance
    >>> covmap, sseg, wseg, fitres = covariance(
    ...     fcube  = 'cube.FILTER_out.fits',
    ...     fvar   = 'varcube.FILTER_out.fits',
    ...     outdir = 'nsigma/',
    ...     mask_source = 'cube.LABELS_out.fits',
    ...     nsegments   = 500,
    ...     plot        = True,
    ... )
    """

    # ------------------------------------------------------------------
    # Set defaults
    # ------------------------------------------------------------------
    if allsizes is None:
        allsizes = list(range(2, 31))

    os.makedirs(outdir, exist_ok=True)

    # ------------------------------------------------------------------
    # Fast path: reload previous sampling results
    # ------------------------------------------------------------------
    if not sampling:
        npzpath = os.path.join(
            outdir, f'npixscaling_spscale_{nsegments}seg.npz'
        )
        if not os.path.isfile(npzpath):
            raise FileNotFoundError(
                f'sampling=False but no previous results found at '
                f'{npzpath}. Run with sampling=True first.'
            )
        if verbose:
            print(f'sampling=False — loading previous results from {npzpath}')
        prev = np.load(npzpath)
        savenormstd = prev['nstd']
        savestd     = prev['std']
        allsizes    = prev['size'].tolist()
        wsegments   = prev['waves']
        nsegments   = len(wsegments) - 1

    # ------------------------------------------------------------------
    # Full sampling path
    # ------------------------------------------------------------------
    if sampling:

        if random_seed is not None:
            np.random.seed(random_seed)
            if verbose:
                print(f'Random seed set to {random_seed}')

        # ------------------------------------------------------------------
        # Open cubes
        # ------------------------------------------------------------------
        if verbose:
            print(f'Opening data cube: {fcube_for_extraction}')
        hducube = fits.open(fcube_for_extraction, memmap=False)
        cube    = hducube[extcube].data
        header  = hducube[extcube].header

        if verbose:
            print(f'Opening variance cube: {fvar_for_extraction}')
        hduvar  = fits.open(fvar_for_extraction, memmap=False)
        vardata = hduvar[extvar].data

        nw, ny, nx = cube.shape
        if verbose:
            print(f'Cube shape: nz={nw}, ny={ny}, nx={nx}')

        # Build wavelength array (support both CD3_3 and CDELT3)
        try:
            delta_lambda = header['CD3_3']
        except KeyError:
            delta_lambda = header['CDELT3']
        zero_lambda = header['CRVAL3']
        wave = np.arange(nw) * delta_lambda + zero_lambda

        # ------------------------------------------------------------------
        # Build the 2-D bad-pixel mask
        # ------------------------------------------------------------------
        # Step 1: edge mask from white-light image
        white_image = np.nanmedian(cube, axis=0)

        if mask_edge is None:
            if verbose:
                print('Building edge mask from white-light image')
            edges   = np.isfinite(white_image)
            badmask = np.ones((ny, nx))
            badmask[edges] = 0.0
            badmask = gaussian_filter(badmask, sigma=1.5)
            badmask[badmask > 0] = 1.0
        else:
            if verbose:
                print(f'Loading edge mask from {mask_edge}')
            badmask = fits.open(mask_edge)[0].data.astype(float)

        # Step 2: source mask
        if mask_source is None:
            if not _SEP_AVAILABLE:
                warnings.warn(
                    'mask_source is None but the sep package is not installed. '
                    'Only the edge mask will be applied. Install sep with '
                    '"pip install sep" to enable automatic source detection.',
                    UserWarning,
                )
            else:
                if verbose:
                    print('Running automatic source detection with sep')
                image_c   = np.ascontiguousarray(white_image, dtype=np.float64)
                badmask_c = np.ascontiguousarray(badmask > 0, dtype=np.uint8)
                bkg       = sep.Background(image_c, mask=badmask_c)
                thresh    = 1.5 * bkg.globalrms
                segmap    = np.zeros((ny, nx))
                _, segmap = sep.extract(
                    image_c - bkg.back(), thresh,
                    segmentation_map=True, minarea=10,
                    clean=True, mask=badmask_c,
                )
                badmask[segmap > 0] = 1.0
        else:
            if verbose:
                print(f'Loading source mask from {mask_source}')
            segmap = fits.open(mask_source)[0].data
            badmask[segmap > 0] = 1.0

        hducube.close()
        hduvar.close()

        # Expand bad-pixel mask to cube shape
        maskcube = np.tile(badmask, (nw, 1, 1))

        # ------------------------------------------------------------------
        # Wavelength segments
        # ------------------------------------------------------------------
        wsegments    = np.linspace(wave[0], wave[-1], nsegments + 1)
        savenormstd  = np.zeros((nsegments, len(allsizes)))
        savestd      = np.zeros((nsegments, len(allsizes)))

        # ------------------------------------------------------------------
        # Helper: build 2-D summed-area tables for a stack of 2-D layers
        # ------------------------------------------------------------------
        def _build_sat(layers):
            """Return the summed-area table for each layer in *layers* (nz, ny, nx).

            The SAT has shape (nz, ny+1, nx+1) so that the sum over the
            rectangle [y0:y1, x0:x1] equals
                sat[:, y1, x1] - sat[:, y0, x1] - sat[:, y1, x0] + sat[:, y0, x0]
            """
            nz, _ny, _nx = layers.shape
            sat = np.zeros((nz, _ny + 1, _nx + 1), dtype=layers.dtype)
            sat[:, 1:, 1:] = np.cumsum(np.cumsum(layers, axis=1), axis=2)
            return sat

        def _sat_rect_sum(sat, y0, y1, x0, x1, ci, dl_loc):
            """Vectorised rectangle sums via summed-area tables.

            Parameters are 1-D int arrays of length N (one per random sample).
            Returns a 1-D float array of length N with the aperture sums.
            """
            # Sum the dl_loc spectral layers using the pre-computed SAT
            total = np.zeros(y0.shape[0], dtype=sat.dtype)
            for k in range(dl_loc):
                layer = ci + k
                total += (  sat[layer, y1, x1]
                          - sat[layer, y0, x1]
                          - sat[layer, y1, x0]
                          + sat[layer, y0, x0])
            return total

        # ------------------------------------------------------------------
        # Main sampling loop
        # ------------------------------------------------------------------
        if verbose:
            print(f'Running covariance sampling loop '
                  f'({nsegments} segments × {len(allsizes)} sizes) …')

        for ii in range(nsegments):
            wmin = wsegments[ii]
            wmax = wsegments[ii + 1]

            index    = np.where((wave > wmin) & (wave <= wmax))[0]
            if len(index) == 0:
                continue
            wavesize = index.max() - index.min()

            pixflux  = cube[index, :, :]
            pixvar   = vardata[index, :, :]
            pixmask  = maskcube[index, :, :]

            if np.all(np.isnan(pixflux)):
                continue

            # Replace NaNs with 0 for SAT computation (NaN would propagate)
            pixflux_safe = np.where(np.isnan(pixflux), 0.0, pixflux)
            pixvar_safe  = np.where(np.isnan(pixvar),  0.0, pixvar)

            # Build summed-area tables once per segment
            sat_flux = _build_sat(pixflux_safe)
            sat_var  = _build_sat(pixvar_safe)
            sat_mask = _build_sat(pixmask)

            # Also build a SAT for NaN counting so we can reject apertures
            # that included any NaN pixel (original code: sum of 0-filled
            # flux and var would silently produce wrong values for such
            # apertures — but the original var would be 0 or negative there,
            # so the normval would be NaN and get filtered out.  We replicate
            # this by tracking where var was originally non-positive.)
            pixvar_bad = np.where((pixvar <= 0) | np.isnan(pixvar), 1.0, 0.0)
            sat_varbad = _build_sat(pixvar_bad)

            ci_max = max(1, wavesize - dl)
            nw_seg = pixflux.shape[0]

            normallstd_seg = []
            allstd_seg     = []

            for size in allsizes:
                # Effective bounds for valid aperture placement
                max_cx = nx - size  # cx can be 0..max_cx (inclusive would give cx+size <= nx)
                max_cy = ny - size

                if max_cx <= 0 or max_cy <= 0:
                    normallstd_seg.append(0.0)
                    allstd_seg.append(0.0)
                    continue

                # --- Batch random draws ---
                all_cx = np.random.uniform(0, nx, size=max_attempts).astype(np.intp)
                all_cy = np.random.uniform(0, ny, size=max_attempts).astype(np.intp)
                all_ci = np.random.uniform(0, ci_max, size=max_attempts).astype(np.intp)

                # Bounds check (vectorised)
                valid = (all_cy + size <= ny) & (all_cx + size <= nx)
                # Also clip ci so ci + dl <= nw_seg
                valid &= (all_ci + dl <= nw_seg)

                idx_valid = np.where(valid)[0]
                if len(idx_valid) == 0:
                    normallstd_seg.append(0.0)
                    allstd_seg.append(0.0)
                    continue

                cx_v = all_cx[idx_valid]
                cy_v = all_cy[idx_valid]
                ci_v = all_ci[idx_valid]

                # SAT rectangle coordinates
                x0 = cx_v
                x1 = cx_v + size
                y0 = cy_v
                y1 = cy_v + size

                # Vectorised aperture sums via SAT
                flux_arr  = _sat_rect_sum(sat_flux, y0, y1, x0, x1, ci_v, dl)
                var_arr   = _sat_rect_sum(sat_var,  y0, y1, x0, x1, ci_v, dl)
                mask_arr  = _sat_rect_sum(sat_mask, y0, y1, x0, x1, ci_v, dl)
                vbad_arr  = _sat_rect_sum(sat_varbad, y0, y1, x0, x1, ci_v, dl)

                # Apply same acceptance criteria as the original code
                accept_mask = (mask_arr < 1) & (var_arr > 0) & (vbad_arr < 1)
                norm_vals   = np.full(len(flux_arr), np.nan)
                ok = accept_mask
                norm_vals[ok] = flux_arr[ok] / np.sqrt(var_arr[ok])

                # Keep only finite values
                finite = np.isfinite(norm_vals)
                good   = np.where(finite)[0]

                # Take at most nsamples (first nsamples that pass, matching
                # the original sequential-draw semantics)
                if len(good) > nsamples:
                    good = good[:nsamples]

                if len(good) > 0:
                    _, _, sigma_norm = sigma_clipped_stats(norm_vals[good])
                    _, _, sigma_flux = sigma_clipped_stats(flux_arr[good])
                else:
                    sigma_norm = 0.0
                    sigma_flux = 0.0

                normallstd_seg.append(sigma_norm)
                allstd_seg.append(sigma_flux)

            savenormstd[ii, :] = normallstd_seg
            savestd[ii, :]     = allstd_seg

            if verbose and (ii % max(1, nsegments // 10) == 0):
                print(f'  Segment {ii + 1}/{nsegments} '
                      f'(λ = {wmin:.1f}–{wmax:.1f} Å)')

        # ------------------------------------------------------------------
        # Save .npz
        # ------------------------------------------------------------------
        npzpath = os.path.join(
            outdir, f'npixscaling_spscale_{nsegments}seg.npz'
        )
        np.savez(
            npzpath,
            nstd=savenormstd,
            std=savestd,
            size=np.array(allsizes),
            waves=wsegments,
        )
        if verbose:
            print(f'Saved covariance data to {npzpath}')

    # ------------------------------------------------------------------
    # Aperture sizes in arcseconds (for plots and fit)
    # ------------------------------------------------------------------
    ssegments = np.array(allsizes) * pixel_scale
    wavecen   = (wsegments[1:] + wsegments[:-1]) / 2.0
    covmap    = savenormstd

    # ------------------------------------------------------------------
    # 1-D polynomial fit per wavelength segment
    # ------------------------------------------------------------------
    fitres = np.zeros((nsegments, fitdeg + 1))
    for ww in range(nsegments):
        try:
            fitres[ww, :] = np.polyfit(ssegments, covmap[ww], fitdeg)
        except Exception:
            pass  # leave zeros if fit fails (e.g., all-zero segment)

    # Build header string for the output text file: Wmin, Wmax, Pn…P0
    poly_header = 'Wmin Wmax ' + ' '.join(
        [f'P{fitdeg - k}' for k in range(fitdeg + 1)]
    )
    output_arr = np.column_stack(
        (wsegments[:-1], wsegments[1:], fitres)
    )
    fitpath = os.path.join(outdir, 'covariance_1Dfit.txt')
    np.savetxt(fitpath, output_arr, header=poly_header)
    if verbose:
        print(f'Saved 1-D polynomial fit to {fitpath}')

    # ------------------------------------------------------------------
    # Diagnostic plot
    # ------------------------------------------------------------------
    if plot:
        matplotlib.rcParams.update({'font.size': 14})

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # --- Left panel: 2-D covariance map ---
        ax = axes[0]
        okcov = (covmap.mean(axis=1) > 0)
        im = ax.pcolor(
            ssegments, wavecen[okcov], covmap[okcov, :],
            vmin=1, vmax=20, shading='nearest', cmap='viridis',
        )
        ax.set_xlabel('Aperture Size (arcsec)')
        ax.set_ylabel('Wavelength (Å)')
        ax.set_title('Empirical covariance map')
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(r'$\sigma_\mathrm{eff}/\sigma_N$')

        # --- Right panel: mean covariance vs aperture + polynomial fit ---
        ax2 = axes[1]
        mean_cov   = np.nanmean(covmap[okcov, :], axis=0)
        fit_coeffs = np.polyfit(ssegments, mean_cov, fitdeg)
        x_fit      = np.linspace(ssegments.min(), ssegments.max(), 200)

        ax2.plot(ssegments, mean_cov, 'ko', label='Data (mean over λ)')
        ax2.plot(
            x_fit, np.polyval(fit_coeffs, x_fit),
            'r-', lw=2,
            label=f'Poly fit (deg={fitdeg})',
        )
        ax2.set_xlabel('Aperture Size (arcsec)')
        ax2.set_ylabel(r'$\sigma_\mathrm{eff}/\sigma_N$')
        ax2.set_title('Mean covariance vs aperture size')
        ax2.legend()

        plt.tight_layout()
        plotpath = os.path.join(outdir, 'covariance_map.pdf')
        plt.savefig(plotpath, bbox_inches='tight')
        plt.close(fig)
        if verbose:
            print(f'Saved covariance plot to {plotpath}')

    return covmap, ssegments, wsegments, fitres
