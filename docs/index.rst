.. SHINE documentation master file.
   Authors: Matteo Fossati, Davide Tornotti

==========================
Documentation of SHINE
==========================

**Authors:** Matteo Fossati, Davide Tornotti

**Date:** 02/07/2026 (v2.1)

.. warning::
   **Version in Development (Unstable Branch)**

   This documentation refers to a version of the code currently under development, restructuring and testing (non-main branch). Some features may not be stable or verified yet.

   For the stable, tested and production-ready version, please refer to the ``main`` branch and the official documentation at `this link <https://shinespec.readthedocs.io/en/latest/>`_.

.. contents::
   :depth: 2
   :local:


Introduction
============

**Project Name:** ``SHINE``

``SHINE`` (Spectral Highlighting and Identification of Emission) is a Python package that identifies connected structures above a user-given signal-to-noise (S/N) threshold in 2D and 3D datasets, and provides a full suite of data-analysis utilities for spectroscopic cubes.

Version 2.1 introduces a modular sub-package architecture:

* ``SHINE`` — Core extraction engine for 2D and 3D data (filtering,
  thresholding, connected-component labelling, catalogue generation).
* ``Find_Em_SHINE`` — End-to-end pipeline for line-emitter detection and
  cataloguing, organised in three steps: extraction, covariance estimation,
  and catalogue building.
* ``Make_Im_SHINE`` — Tool to generate 2-D surface-brightness images from
  spectroscopic cubes using the SHINE extraction products.

Previous changes (v2.0):

* ``shine_utils`` — general-purpose utility functions (continuum subtraction).
* Removal of the Tkinter GUI (``GUI_SHINE``).
* ``matplotlib`` added as a core dependency; ``sep`` added as an optional
  dependency for automatic source masking in the covariance estimator.
* Fully revised NumPy-style docstrings for all public functions.


Installation
============

**Requirements:**

- Python version: ``>=3.8``
- Core dependencies: ``numpy``, ``scipy``, ``astropy``, ``matplotlib``, ``connected-components-3d``
- Optional: ``sep`` (for automatic continuum-source detection in :func:`~shine.Find_Em_SHINE.covariance`)
- Optional: ``mypython`` (for 1-D spectral extraction in :func:`~shine.Find_Em_SHINE.build_em_catalog`)

**Steps to Install:**

1. Clone the repository::

       git clone https://github.com/matteofox/SHINE.git
       cd SHINE

2. Install the code::

       python -m pip install .

3. (Optional) Install with covariance support::

       python -m pip install ".[covariance]"


Package Structure
==================

The ``shine/`` package is organised as follows::

    shine/
    ├── __init__.py             # Package root: exposes sub-packages and runextraction
    ├── SHINE.py                # Core extraction engine
    ├── shine_utils.py          # General-purpose utilities (clean_clube)
    │
    ├── Find_Em_SHINE/          # Line-emitter extraction and cataloguing
    │   ├── __init__.py         # Exposes: extract, covariance, build_em_catalog
    │   ├── extraction.py       # Step 1: SHINE extraction with emitter defaults
    │   ├── covariance.py       # Step 2: empirical noise-covariance estimation
    │   └── catalogue.py        # Step 3: final catalogue, S/N correction, cutouts
    │
    └── Make_Im_SHINE/          # 2-D image generation
        ├── __init__.py         # Exposes: Make_Im_SHINE, main
        └── Make_Im_SHINE.py    # Image creation from cubes and labels


Usage of SHINE — Core Extraction
=================================

The extraction is performed using ``SHINE``. The basic idea behind the code is as follows:

1. (Optional, applicable to 3D data only) Select a portion of the cube (in the z-direction or wavelength direction) where the user wants to focus the extraction.
2. Mask certain voxels using a user-provided mask (e.g., continuum sources).
3. Spatially filter the cube/image and the associated 2D or 3D variance using a user-defined kernel dimension.
4. Apply a threshold to the cube/image based on the user-defined S/N threshold.
5. Group connected voxels (3D)/pixels (2D) that meet the S/N threshold and other user-defined parameters.
6. Generate and save the catalog along with the labeled cube/image.


Run Extraction from the command line
--------------------------------------

**Basic Usage for 3D cube:**

.. code-block:: bash

   SHINE <3D-cube> <3D-variance>

- **Example (without selecting a subcube):**

.. code-block:: bash

   SHINE ../Cubes/Datacube.fits ../Cubes/Datavarcube.fits --snthreshold 2 --spatsmooth 3 --minvox 300 --mindz 2 --outdir ../Dataproducts/ --writelabels

- **Example (selecting a subcube):**

.. code-block:: bash

   SHINE ../Cubes/Datacube.fits ../Cubes/Datavarcube.fits --zmin 40 --zmax 100 --snthreshold 2 --spatsmooth 3 --minvox 300 --mindz 2  --outdir ../Dataproducts/ --writelabels --writesubcube


**Basic Usage for 2D images:**

.. code-block:: bash

   SHINE <2D-image> <2D-variance>

- **Example:**

.. code-block:: bash

   SHINE ../Cubes/Dataimage.fits ../Cubes/Datavarimage.fits --snthreshold 2 --spatsmooth 3 --minvox 300 --outdir ../Dataproducts/ --writelabels



**Command line arguments for SHINE:**


*General Arguments:*

- ``-h, --help``: Show this help message and exit.

*Input Control Arguments:*

- ``data``: Path of the input data (3D or 2D). Expected to be in extension 0, unless ``extdata`` is defined.
- ``vardata``: Path of the variance cube. Expected to be in extension 0, unless ``extvardata`` is defined. If a positive number is provided a constant variance is assumed; if -1, the variance is computed from the data with a sigma-clipping algorithm.
- ``--mask2d``: (Optional) Path of an optional two-dimensional mask to be applied along the wave axis.
- ``--mask2dpost``: (Optional) Path of an optional two-dimensional mask to be applied after the spatial smoothing.
- ``--mask3d``: (Optional) Path of an optional three-dimensional mask. **(Not implemented yet)**.
- ``--extdata``: Specifies the HDU index in the FITS file to use for data extraction (default=0).
- ``--extvardata``: Specifies the HDU index in the FITS file variance to use for data extraction (default=0).
- ``--zmin``: (Optional) Select the cube and the variance: initial pixel in z direction (from 0). Only valid for 3D data.
- ``--zmax``: (Optional) Select the cube and the variance: final pixel in z direction (from 0). Only valid for 3D data.
- ``--lmin``: (Optional) Select the cube and the variance: initial wavelength in z direction (in Angstrom). Only valid for 3D data.
- ``--lmax``: (Optional) Select the cube and the variance: final wavelength in z direction (in Angstrom). Only valid for 3D data.

*Extraction Arguments:*

- ``--snthreshold``: The SNR of voxels (3D)/pixels (2D) to be included in the extraction (default=2).
- ``--spatsmooth``: Gaussian Sigma of the spatial convolution kernel applied in X and Y (default=0).
- ``--spatsmoothX``: (Optional) Gaussian Sigma of the spatial convolution kernel applied in X. Has priority over ``spatsmooth``.
- ``--spatsmoothY``: (Optional) Gaussian Sigma of the spatial convolution kernel applied in Y. Has priority over ``spatsmooth``.
- ``--specsmooth``: Gaussian Sigma of the spectral convolution kernel applied in Lambda.
- ``--usefftconv``: If ``True``, use FFT for convolution rather than the direct algorithm.
- ``--dovarsmooth``: If False, do not apply the smoothing on the vardata.
- ``--connectivity``: Voxel connectivity scheme to be used (default=26). Allowed values: 4, 8 (2D); 26, 18, 6 (3D).
- ``--maskspedge``: Determines how much, in pixels (default=20), to expand the mask around the edges of the cube/image.

*Cleaning Arguments:*

- ``--minvox``: Minimum number of connected voxels (3D)/pixels (2D) for a source to be in the final catalogue (default=1). For 2D data this argument has priority over ``--minarea``.
- ``--mindz``: Minimum number of connected voxels in the spectral direction for a source to be in the final catalogue (default=1). Only valid for 3D data.
- ``--maxdz``: Maximum number of connected voxels in the spectral direction for a source to be in the final catalogue (default=200). Only valid for 3D data.
- ``--minarea``: Minimum number of connected projected spatial voxels (3D)/pixels (2D) for a source to be in the final catalogue (default=1).

*Output Control Arguments (inherit the original file name):*

- ``--outdir``: Output directory path (Default ./).
- ``--writelabels``: If set, write labels cube/image. The file is saved as ``dataname.LABELS_out.fits``.
- ``--writesmdata``: If set, write the smoothed cube/image. The file is saved as ``dataname.FILTER_out.fits``.
- ``--writesmvar``: If set, write the smoothed variance. The file is saved as ``varname.FILTER_out.fits``.
- ``--writesmsnr``: If set, write the S/N smoothed cube/image. The file is saved as ``dataname.FILTERSNR_out.fits``.
- ``--writesubcube``: If set and used, write the subcubes (cube and variance). Only valid for 3D data. The file is saved as ``dataname.SUBCUBE.fits`` or ``varname.SUBCUBE.fits``.
- ``--writevardata``: If vardata is a user-provided constant value or -1 (auto-computed from sigma-clipping) write the variance.


Run Extraction using Python
-----------------------------

SHINE can be used also in a Python code.

**Basic Usage for 3D data:**

.. code-block:: python

    from astropy.table import Table
    from shine import SHINE

    # Extraction on the full cube
    SHINE.runextraction(
        '../Data/Datacube.fits', '../Data/Datavarcube.fits',
        snthreshold=2, spatsmooth=4,
        minvox=3000, minarea=1000,
        mask2d='../Data/2D_MASK.fits',
        mask2dpost='../Data/2D_MASK_post.fits',
        outdir='../Dataproducts/',
        writelabels=True, maskspedge=20,
        writesmdata=True, writesmsnr=True,
    )

    # Extraction on a subcube (by pixel index)
    SHINE.runextraction(
        '../Data/Datacube.fits', '../Data/Datavarcube.fits',
        zmin=40, zmax=100,
        snthreshold=2, spatsmooth=4,
        minvox=3000, minarea=1000, maskspedge=20,
        mask2d='../Data/2D_MASK.fits',
        mask2dpost='../Data/2D_MASK_post.fits',
        outdir='../Dataproducts/',
        writelabels=True, writesmdata=True, writesmsnr=True,
    )

    # Quick visualisation of the output catalogue
    catalogue3D = Table.read('../Outdir/Datacube.CATALOGUE_out.fits')



.. figure:: /_static/Catalogue_3Dextraction.png
   :width: 80%
   :align: left
   :alt: Catalogue output from 3D extraction.

   Catalogue output from 3D extraction (Example).

**Basic Usage for 2D data:**

.. code-block:: python

    from astropy.table import Table
    from shine import SHINE

    # Attention! Remember to change connectivity to 4 or 8 for 2D data.
    SHINE.runextraction(
        '../Data/Dataimage.fits', '../Data/Datavarimage.fits',
        connectivity=8, snthreshold=3, spatsmooth=1,
        minvox=40, maskspedge=20,
        outdir='../Dataproducts/', writelabels=True,
    )

    catalogue2D = Table.read('../Outdir/Dataimage.CATALOGUE_out.fits')



.. figure:: /_static/Catalogue_2Dextraction.png
   :width: 80%
   :align: left
   :alt: Catalogue output from 2D extraction.

   Catalogue output from 2D extraction (Example).


Usage of Find_Em_SHINE — Line-Emitter Pipeline
================================================

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

   # products is a dictionary with paths to the output files:
   # products['fcube_filtered']  — smoothed data cube
   # products['fvar_filtered']   — smoothed variance cube
   # products['fsegmap']         — segmentation (labels) cube
   # products['fcatalogue']      — raw extraction catalogue
   # products['fcube_clean']     — continuum-subtracted cube (if do_continuum_sub=True)


**Parameters:**

- ``fcube``: Path to the science data cube FITS file.
- ``fvar``: Path to the variance cube FITS file, or a numeric string / ``'-1'``.
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
       fcube       = products['fcube_filtered'],
       fvar        = products['fvar_filtered'],
       outdir      = './covariance/',
       mask_source = products['fsegmap'],     # optional
       nsegments   = 500,
       allsizes    = list(range(2, 31)),
       dl          = 4,
       nsamples    = 10000,
       pixel_scale = 0.2,
       fitdeg      = 2,
       plot        = True,
   )

**Parameters:**

- ``fcube``: Path to the filtered data cube (``*FILTER_out.fits`` from SHINE ``--writesmdata``).
- ``fvar``: Path to the filtered variance cube (``*FILTER_out.fits`` from SHINE ``--writesmvar``).
- ``outdir``: Output directory (created automatically if absent).
- ``extcube``, ``extvar`` *(default 0)*: HDU extension indices.
- ``allsizes`` *(default range(2, 31))*: List of aperture sizes in pixels to probe.
- ``nsegments`` *(default 500)*: Number of spectral wavelength segments.
- ``dl`` *(default 4)*: Spectral depth (layers) of each random aperture sample.
- ``nsamples`` *(default 10 000)*: Valid samples to collect per (size, segment) combination.
- ``max_attempts`` *(default 100 000)*: Maximum random draws before giving up on a combination.
- ``mask_source`` *(default None)*: Path to a FITS source/continuum mask. If *None* and ``sep`` is installed, sources are detected automatically; otherwise only the edge mask is applied.
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
       fcube       = products['fcube_filtered'],
       fcube_var   = products['fvar_filtered'],
       fsegmap     = products['fsegmap'],
       catpath     = products['fcatalogue'],
       outdir      = './emitters/',
       cov_dir     = './covariance/',  
       SNcut       = (7, 5),
       checkimg    = True,
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
8. If ``fcube_orig`` is provided and ``mypython`` is installed, extracts
   1-D spectra for each source.

**Parameters:**

- ``fcube``: Path to the filtered data cube.
- ``fcube_var``: Path to the filtered variance cube.
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
- ``fcube_orig`` *(default None)*: Unsmoothed cube for 1-D spectral extraction (requires ``mypython``).
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
   import numpy as np

   # ---- Step 1: Extract ----
   products = shine.Find_Em_SHINE.extract(
       'Datacube.fits', 'Varcube.fits',
       snthreshold=2.0, spatsmooth=2.0,
       connectivity=26, maskspedge=20,
       mindz=1, maxdz=200, minvox=1, minarea=1,
       outdir='./extraction/',
       do_continuum_sub=True,
   )

   # ---- Step 2: Covariance ----
   covmap, sseg, wseg, fitres = shine.Find_Em_SHINE.covariance(
       fcube       = products['fcube_filtered'],
       fvar        = products['fvar_filtered'],
       outdir      = './covariance/',
       mask_source = products['fsegmap'],
       nsegments   = 500,
       pixel_scale = 0.2,
   )

   # ---- Step 3: Build catalogue ----
   catalog = shine.Find_Em_SHINE.build_em_catalog(
       fcube     = products['fcube_filtered'],
       fcube_var = products['fvar_filtered'],
       fsegmap   = products['fsegmap'],
       catpath   = products['fcatalogue'],
       outdir    = './emitters/',
       cov_dir   = './covariance/',   # automatically loads covariance_1Dfit.txt
       target_z  = 3.1,
       rest_line = 1216.0,
       vel_cut   = 2000,
       SNcut     = (7, 5),
       checkimg  = True,
       padding   = 50,
   )

   print(f'Final catalogue: {len(catalog)} sources')


Usage of Make_Im_SHINE — Image Generation
==========================================

``Make_Im_SHINE`` creates 2-D surface-brightness images from 3-D spectroscopic
cubes, optionally using the SHINE labels cube to select specific sources.


Run generation of images by command line
-----------------------------------------

**Basic Usage (for 3D data only):**

.. code-block:: bash

   Make_Im_SHINE <Datacube> <Cubelabels>

- **Example using a 3D mask for selecting voxels:**

.. code-block:: bash

   Make_Im_SHINE ../Cubes/Datacube.fits ../Cubes/Cubelabels.fits --Id 2 5 9  --outdir ../Dataproducts/ --itype flux --writeout

- **Example using a 3D mask for selecting one single object and the associated pseudo narrow-band around it (only for one single Id):**

.. code-block:: bash

   Make_Im_SHINE ../Cubes/Datacube.fits ../Cubes/Cubelabels.fits --Id 2  --outdir ../Dataproducts/ --itype flux --nls -1 --nlsadd 2 --writeout

- **Example to create a narrow band image:**

.. code-block:: bash

    Make_Im_SHINE ../Cubes/Datacube.fits  --outdir ../Dataproducts/ --itype flux --writeout


**Command line arguments for Make_Im_SHINE:**

*General Arguments:*

- ``-h, --help``: Show this help message and exit.

*Input Control Arguments:*

- ``cube``: Path of the input datacube. Expected to be in extension 0, unless ``extcub`` is defined.
- ``labelsCube``: Path of the cube with labels. Expected to be in extension 0, unless ``extlabels`` is defined.
- ``--Id``: The IDs of the grouped voxels to be used for the surface brightness image extraction. If a list is passed, all the valid IDs will be stacked into the final image. If ``[-1]`` (default), stacks all the IDs in the ``labelsCube``.
- ``--itype``: Type of image to produce. Allowed values: ``flux`` (default), ``mean``, or ``median``.
- ``--extcub``: Specifies the HDU index in the FITS file cube to use for data extraction (default=0).
- ``--extlabels``: Specifies the HDU index in the FITS file labels to use for labels data extraction (default=0).
- ``--nsl``: Pseudo narrow-band: selects the provided layer associated with the object as the central layer from which the pseudo narrow-band is built. If ``-1``, it selects the mean layer of the object. Use this only if ``len(Id)=1``. Default is ``-2`` (no noise layers).
- ``--nsladd``: Pseudo narrow-band: specifies how many layers to collapse adjacent to the selected central one (default=0).

*Output Control Arguments:*

- ``--outdir``: Output directory path (default=``./``).
- ``--writeout``: If set, writes the flux image and the associated variance image. The file is saved as ``dataname.IMAGE.fits``.
- ``--addname``: Optional suffix to append to the base name of the output file.


Run generation of images using Python
--------------------------------------

.. code-block:: python

   from shine.Make_Im_SHINE import Make_Im_SHINE

   # Image from labels cube, single ID, no pseudo narrow-band
   img = Make_Im_SHINE(
       '../Data/Datacube.fits', labelsCube='../Data/Labelscube.fits',
       Id=[45], extcub=0, extlabels=0, itype='flux',
       outdir='../Dataproducts', writeout=True,
   )

   # Image from labels cube, single ID, pseudo narrow-band ±2 layers around mean
   img = Make_Im_SHINE(
       '../Data/Datacube.fits', labelsCube='../Data/Labelscube.fits',
       Id=[45], extcub=0, extlabels=0, itype='flux',
       outdir='../Dataproducts', writeout=True, nsl=-1, nsladd=2,
   )

   # Narrow-band image (all voxels)
   img = Make_Im_SHINE(
       '../Data/Datacube.fits', extcub=0, itype='flux',
       outdir='../Dataproducts', writeout=True,
   )


General Utilities (``shine_utils``)
=====================================

The ``shine_utils`` module provides general-purpose utility functions.


clean_clube — continuum subtraction
-------------------------------------

.. code-block:: python

   from shine.shine_utils import clean_clube
   from astropy.io import fits

   cube = fits.open('Datacube.fits')[0].data
   cube_clean = clean_clube(cube, filtsize=7, rebinfac=40)

**Parameters:**

- ``data``: 3-D numpy array (nz, ny, nx).
- ``filtsize`` *(default 7)*: Number of re-binned spectral slices used for the median filter applied along the spectral axis.
- ``rebinfac`` *(default 40)*: Number of original spectral layers collapsed into each continuum slice.

**Returns:** continuum-subtracted masked array with the same shape as the input.

**Algorithm:**

1. Re-bin the cube spectrally into groups of ``rebinfac`` layers using a sigma-clipped median.
2. Apply a 1-D median filter of width ``filtsize`` along the re-binned spectral axis.
3. Subtract the interpolated continuum layer-by-layer from the original cube.

.. note::
   ``clean_clube`` can also be called automatically as part of
   ``Find_Em_SHINE.extract()`` by setting ``do_continuum_sub=True``.


filter_cube — spatial and spectral smoothing
----------------------------------------------

.. code-block:: python

   from shine.shine_utils import filter_cube

   # Spatial smoothing on a 3-D cube
   cube_smoothed = filter_cube(cube, spatsmooth=2.0, specsig=0.0)

Applies a 2-D or 3-D Gaussian smoothing kernel to a data cube or variance image.

**Parameters:**

- ``cube``: Input 2-D image or 3-D data cube.
- ``spatsmooth`` *(default 2.0)*: Gaussian sigma (pixels) for spatial smoothing. Can be a single float or list of two floats ``[sig_x, sig_y]``.
- ``specsig`` *(default 0.0)*: Gaussian sigma (pixels) for spectral axis smoothing. Only valid for 3-D data.
- ``isvar`` *(default False)*: If ``True``, the kernel is squared and not normalized to unity (tailored for variance data).
- ``usefftconv`` *(default False)*: Use FFT-based convolution instead of direct convolution.


subcube — spectral sub-cube extraction
----------------------------------------

.. code-block:: python

   from shine.shine_utils import subcube

   # Extract a sub-cube by pixel layers
   sub_cube, sub_header = subcube(cube, datahead=header, zmin=40, zmax=100)

Extracts a spectral sub-cube by pixel index (Z axis) or wavelength range.

**Parameters:**

- ``cube``: Input 3-D data cube array.
- ``datahead``: FITS header of the input cube (used to extract wavelength WCS keywords).
- ``filename``: Stem name of the file (used for writing the output file).
- ``pathcube``: Path to the FITS file of the cube. If provided, `cube`, `datahead`, and `filename` are read directly from it.
- ``extcube`` *(default 0)*: HDU extension index to read from `pathcube`.
- ``zmin``, ``zmax``: Spectral layer pixel index range (0-based) to extract.
- ``lmin``, ``lmax``: Wavelength range (in Å) to extract.
- ``writesubcube`` *(default False)*: Write the sub-cube to a FITS file.
- ``addname``: Suffix to append to the output filename.


.. _changelog:

Changelog
=========

.. include:: ../CHANGELOG


Contributing
============

If you are interested in contributing to the project, please contact us and follow these steps:

1. Fork the repository on GitHub.
2. Create a new branch for your feature/bugfix.
3. Submit a pull request.


API
===

This section lists the main available functions in the SHINE package.

**shine.SHINE** — Core Extraction Engine

- ``SHINE.runextraction``: Performs the full extraction process (filter → threshold → CC labelling → catalogue → photometry → WCS).
- ``SHINE.threshold_cube``: Thresholds the cube in S/N, returning a boolean cube.
- ``SHINE.extract``: Applies the CC3D connected-components algorithm.
- ``SHINE.generate_catalogue``: Generates a catalogue from a labelled cube.
- ``SHINE.cleaning``: Removes sources that do not meet size/voxel criteria.
- ``SHINE.compute_photometry``: Measures flux, flux uncertainty, S/N, and flux-weighted centroids.
- ``SHINE.add_wcs_struct``: Adds RA/Dec/Lambda columns to the catalogue.

**shine.Find_Em_SHINE** — Line-Emitter Pipeline

- ``Find_Em_SHINE.extract``: Run SHINE extraction with emitter-optimised defaults; optional continuum subtraction.
- ``Find_Em_SHINE.covariance``: Empirical noise-covariance estimation with random aperture sampling.
- ``Find_Em_SHINE.build_em_catalog``: Full post-processing pipeline: S/N correction, confidence classes, cutouts, spectra.

**shine.Make_Im_SHINE** — Image Generation

- ``Make_Im_SHINE.Make_Im_SHINE``: Produces 2-D surface-brightness images (flux/mean/median) from the cube using the labels cube.

**shine.shine_utils** — General Utilities

- ``shine_utils.clean_clube``: Continuum subtraction via spectral re-binning and median filtering.
- ``shine_utils.filter_cube``: Spatial and spectral Gaussian smoothing.
- ``shine_utils.subcube``: Spectral sub-cube extraction.


License
=======

Copyright (C) 2024–2026 The Authors

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License.
