Line-Emitter Pipeline (``Find_Em_SHINE``)
=========================================

``Find_Em_SHINE`` provides a three-step modular pipeline for extracting and
cataloguing line emitters from spectroscopic cubes. Each step can be called
independently, enabling fine-grained control over the workflow.

.. code-block:: python

   import shine

   # Step 1: Extraction
   products = shine.Find_Em_SHINE.extract(...)

   # Step 2: Covariance estimation
   covmap, sseg, wseg, fitres = shine.Find_Em_SHINE.covariance(...)

   # Step 3: Catalogue building
   catalog = shine.Find_Em_SHINE.build_em_catalog(...)


Run Emitter Pipeline from the command line
--------------------------------------------

You can run the entire 3-step pipeline (extraction, covariance, and catalogue) with a single command from your terminal:

.. code-block:: bash

   Find_Em_SHINE <path/to/config.ini>

**Initialize a template configuration file:**

If you don't have a configuration file yet, you can generate a default template containing all parameters and helpful comments:

.. code-block:: bash

   Find_Em_SHINE --init-config my_config.ini

Open the generated file, edit the paths and parameters under the sections ``[FILES]``, ``[CONTINUUM_SUBTRACTION]``, ``[EXTRACTION]``, ``[COVARIANCE]``, and ``[CATALOGUE]``, and run the pipeline.


Step 1 — ``extract``: SHINE extraction for emitters
------------------------------------------------------

A convenience wrapper around ``SHINE.runextraction()`` with defaults optimised
for line-emitter detection. It:

* optionally performs continuum subtraction before extraction
  (via ``clean_clube``);
* forces output of labels, filtered data, and filtered variance cubes
  (required by the downstream covariance and catalogue steps).

.. code-block:: python

   from shine.Find_Em_SHINE import extract

   products = extract(
       '../Data/Datacube.fits',
       '../Data/Datavarcube.fits',
       snthreshold=3.0,
       spatsmooth=2.0,
       connectivity=26,
       maskspedge=0,
       mindz=3,
       maxdz=50,
       minvox=27,
       minarea=9,
       outdir='./extraction/',
       do_continuum_sub=True,    # optional continuum subtraction
       rebinfac=40,
       filtsize=7,
   )

**Parameters:**

- ``fcube_input``: Path to the science data cube FITS file.
- ``fvar_input``: Path to the variance cube FITS file, or a numeric string / ``'-1'``.
- ``extdata``, ``extvar`` *(default 0)*: HDU extension indices.
- ``mask2d``, ``mask2dpost`` *(default None)*: Paths to pre/post-smoothing 2-D masks (1 = bad).
- ``snthreshold`` *(default 2.0)*: S/N threshold for voxel inclusion.
- ``spatsmooth`` *(default 2.0)*: Spatial Gaussian sigma (pixels).
- ``specsmooth`` *(default 0.0)*: Spectral Gaussian sigma (pixels); 0 to disable.
- ``connectivity`` *(default 26)*: Voxel connectivity: 6, 18, or 26.
- ``maskspedge`` *(default 20)*: Pixels to mask around the field edges.
- ``mindz`` *(default 1)*: Minimum spectral extent per source (layers).
- ``maxdz`` *(default 200)*: Maximum spectral extent per source (layers).
- ``minvox`` *(default 1)*: Minimum total voxels per source.
- ``minarea`` *(default 1)*: Minimum projected spatial area per source.
- ``zmin``, ``zmax`` *(default None)*: Layer index range.
- ``lmin``, ``lmax`` *(default None)*: Wavelength range (Å).
- ``outdir`` *(default './')*: Output directory.
- ``do_continuum_sub`` *(default False)*: If ``True``, subtract the continuum before extraction using ``clean_clube``.
- ``rebinfac`` *(default 40)*: Spectral re-binning factor for continuum subtraction.
- ``filtsize`` *(default 7)*: Median filter width for continuum subtraction.

**Returns:** a dictionary with keys ``'fcube_filtered'``, ``'fvar_filtered'``,
``'fsegmap'``, ``'fcatalogue'``, ``'fcube_clean'``.


Step 2 — ``covariance``: empirical noise-covariance estimation
-----------------------------------------------------------------

Estimates the empirical noise covariance of the (smoothed) data cube by
sampling random apertures of increasing projected size at random positions
and wavelengths. The measured flux/noise ratio is fitted with a per-wavelength
polynomial model that can later be used to correct the S/N of extracted sources.

.. code-block:: python

   from shine.Find_Em_SHINE import covariance

   covmap, sseg, wseg, fitres = covariance(
       fcube_for_extraction = products['fcube_filtered'],
       fvar_for_extraction  = products['fvar_filtered'],
       outdir               = './covariance/',
       nsegments            = 500,
       allsizes             = list(range(2, 31)),
       dl                   = 4,
       nsamples             = 10000,
       pixel_scale          = 0.2,
       fitdeg               = 2,
       plot                 = True,
   )

**Parameters:**

- ``fcube_for_extraction``: Path to the filtered data cube (``*FILTER_out.fits`` from SHINE ``--writesmdata``).
- ``fvar_for_extraction``: Path to the filtered variance cube (``*FILTER_out.fits`` from SHINE ``--writesmvar``).
- ``outdir``: Output directory (created automatically if absent).
- ``extcube``, ``extvar`` *(default 0)*: HDU extension indices.
- ``allsizes`` *(default range(2, 31))*: List of aperture sizes in pixels to probe.
- ``nsegments`` *(default 500)*: Number of spectral wavelength segments.
- ``dl`` *(default 4)*: Spectral depth (layers) of each random aperture sample.
- ``nsamples`` *(default 10 000)*: Valid samples to collect per (size, segment) combination.
- ``max_attempts`` *(default 100 000)*: Maximum random draws before giving up on a combination.
- ``mask_source`` *(default None)*: Path to a 2D FITS source/continuum mask. If *None* and ``sep`` is installed, sources are detected automatically; otherwise only the edge mask is applied.
- ``mask_edge`` *(default None)*: Path to a FITS edge mask. If *None*, derived automatically from NaN/zero pixels.
- ``pixel_scale`` *(default 0.2)*: Pixel scale in arcsec/pixel (MUSE nominal).
- ``fitdeg`` *(default 2)*: Degree of the polynomial fit to the covariance per wavelength segment.
- ``plot`` *(default True)*: Save a diagnostic 2-panel PDF.
- ``verbose`` *(default True)*: Print progress.

**Output files written to** ``outdir``:

- ``npixscaling_spscale_{nsegments}seg.npz`` — NumPy archive with keys ``nstd``, ``std``, ``size``, ``waves``.
- ``covariance_map.pdf`` — 2-panel diagnostic plot (2-D covariance map + mean covariance vs aperture size with polynomial fit).
- ``covariance_1Dfit.txt`` — ASCII table with columns ``Wmin Wmax P2 P1 P0`` (columns scale with ``fitdeg``).

**Returns:** ``(covmap, ssegments, wsegments, fitres)`` tuple.


Step 3 — ``build_em_catalog``: final emitter catalogue
---------------------------------------------------------

Builds the final catalogue of line emitters starting from the raw SHINE
catalogue and segmentation map. Computes covariance-corrected S/N, applies
quality cuts, assigns confidence classes, and generates per-source image
cutouts and spectra.

.. code-block:: python

   from shine.Find_Em_SHINE import build_em_catalog

   catalog = build_em_catalog(
       fcube_for_extraction = products['fcube_filtered'],
       fvar_for_extraction  = products['fvar_filtered'],
       fsegmap              = products['fsegmap'],
       catpath              = products['fcatalogue'],
       outdir               = './emitters/',
       cov_dir              = './covariance/',
       SNcut                = (7, 5),
       checkimg             = True,
   )

**Main steps performed:**

1. Loads the raw SHINE catalogue and adds quality columns.
2. Computes the covariance correction vector from ``cov_poly``.
3. Computes velocity offsets and (optionally) trims the catalogue.
4. Computes per-source corrected S/N, filling ``SNR``, ``SNR_odd``,
   ``SNR_even``, ``SNR_med``, ``covfac``, ``BoxFraction``, ``OverContinuum``.
5. Applies S/N and quality cuts and assigns confidence classes.
6. Writes ``{catname}_all_SNR.fits`` and ``{catname}_select_SNR.fits``.
7. If ``checkimg=True``, generates per-source image cutouts.
8. If ``fcube_for_spectra`` is provided and ``mypython`` is installed, extracts
   1-D spectra for each source.

**Parameters:**

- ``fcube_for_extraction``: Path to the filtered data cube.
- ``fvar_for_extraction``: Path to the filtered variance cube.
- ``fsegmap``: Path to the SHINE segmentation map (``*LABELS_out.fits``).
- ``catpath``: Path to the SHINE output catalogue (``*CATALOGUE_out.fits``).
- ``outdir`` *(default './')*: Output directory.
- ``cov_poly`` *(default None)*: Covariance correction model. Accepts a 1-D array (single polynomial in aperture size) or a 2-D array (per-wavelength polynomial, as written by ``covariance``). If ``None``, no correction is applied.
- ``cov_dir`` *(default None)*: Path to the covariance estimation directory. If provided and ``cov_poly`` is ``None``, the 1-D polynomial fit coefficients are automatically loaded from ``{cov_dir}/covariance_1Dfit.txt``.
- ``target_z``, ``rest_line`` *(default None)*: Target redshift and rest-frame line wavelength for velocity-offset computation.
- ``vel_cut`` *(default None)*: Remove sources with ``|veloffset| > vel_cut`` (km/s).
- ``SNcut`` *(default (7, 5))*: S/N thresholds for confidence class assignment (highest first).
- ``DeltaEOSNcut`` *(default (0.5, 0.5))*: Max fractional even/odd S/N difference per class.
- ``SNEOcut`` *(default (3, 3))*: Minimum even and odd S/N per class.
- ``fracnpix`` *(default None)*: Minimum box fraction for source inclusion.
- ``derived`` *(default True)*: Compute S/N and other derived quantities.
- ``checkimg`` *(default True)*: Generate per-source image cutouts.
- ``mask`` *(default None)*: Path to a 2-D mask (1 = bad) for removing sources.
- ``startind`` *(default 0)*: Index from which to start the cutout loop.
- ``padding`` *(default 50)*: Spatial padding (pixels) for image cutouts.
- ``pixel_scale`` *(default 0.2)*: Pixel scale in arcsec/pixel.
- ``fcube_median``, ``fcube_odd``, ``fcube_even`` *(default None)*: Half-exposure cubes for quality checks.
- ``fcube_median_var``, ``fcube_odd_var``, ``fcube_even_var`` *(default None)*: Associated variance cubes.
- ``fcube_for_spectra`` *(default None)*: Unsmoothed cube for 1-D spectral extraction (requires ``mypython``).
- ``fsource_img`` *(default None)*: 2-D continuum aperture map for overlap flagging.
- ``marzred`` *(default None)*: Marz-format redshift catalogue for continuum sources.

**Output files:**

- ``{outdir}/{catname}_all_SNR.fits`` — all sources above the minimum S/N.
- ``{outdir}/{catname}_select_SNR.fits`` — sources after all quality cuts.
- ``{outdir}/objs/id{Id}/id{Id}_img.fits`` — per-source image cutouts.
- ``{outdir}/objs/id{Id}/spectrum.fits`` — 1-D spectra (requires ``mypython``).

**Returns:** ``astropy.table.Table`` — the final catalogue.


Full Find_Em_SHINE Workflow Example
---------------------------------------

Here is a complete example showing how to run the full emitter extraction
pipeline from Python:

.. code-block:: python

   import shine

   # ---- Step 1: Extract ----
   products = shine.Find_Em_SHINE.extract(
       fcube_input='Datacube.fits', fvar_input='Varcube.fits',
       snthreshold=2.0, spatsmooth=2.0,
       connectivity=26, maskspedge=20,
       mindz=1, maxdz=200, minvox=1, minarea=1,
       outdir='./extraction/',
       do_continuum_sub=True,
   )

   # ---- Step 2: Covariance ----
   covmap, sseg, wseg, fitres = shine.Find_Em_SHINE.covariance(
       fcube_for_extraction = products['fcube_filtered'],
       fvar_for_extraction  = products['fvar_filtered'],
       outdir               = './covariance/',
       mask_source          = products['fsegmap'],
       nsegments            = 500,
       pixel_scale          = 0.2,
   )

   # ---- Step 3: Build catalogue ----
   catalog = shine.Find_Em_SHINE.build_em_catalog(
       fcube_for_extraction = products['fcube_filtered'],
       fvar_for_extraction  = products['fvar_filtered'],
       fsegmap              = products['fsegmap'],
       catpath              = products['fcatalogue'],
       outdir               = './emitters/',
       cov_dir              = './covariance/',   # automatically loads covariance_1Dfit.txt
       target_z             = 3.1,
       rest_line            = 1216.0,
       vel_cut              = 2000,
       SNcut                = (7, 5),
       checkimg             = True,
       padding              = 50,
   )

   print(f'Final catalogue: {len(catalog)} sources')
