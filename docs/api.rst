API Reference
=============

This section lists the main available modules and functions in the SHINE package.

**shine.SHINE** — Core Extraction Engine
----------------------------------------

- ``SHINE.runextraction``: Performs the full extraction process (filter → threshold → CC labelling → catalogue → photometry → WCS).
- ``SHINE.threshold_cube``: Thresholds the cube in S/N, returning a boolean cube.
- ``SHINE.extract``: Applies the CC3D connected-components algorithm.
- ``SHINE.generate_catalogue``: Generates a catalogue from a labelled cube.
- ``SHINE.cleaning``: Removes sources that do not meet size/voxel criteria.
- ``SHINE.compute_photometry``: Measures flux, flux uncertainty, S/N, and flux-weighted centroids.
- ``SHINE.add_wcs_struct``: Adds RA/Dec/Lambda columns to the catalogue.

**shine.Find_Em_SHINE** — Line-Emitter Pipeline
-----------------------------------------------

- ``Find_Em_SHINE.extract``: Run SHINE extraction with emitter-optimised defaults; optional continuum subtraction.
- ``Find_Em_SHINE.covariance``: Empirical noise-covariance estimation with random aperture sampling.
- ``Find_Em_SHINE.build_em_catalog``: Full post-processing pipeline: S/N correction, confidence classes, cutouts, spectra.

**shine.Make_Im_SHINE** — Image Generation
------------------------------------------

- ``Make_Im_SHINE.Make_Im_SHINE``: Produces 2-D surface-brightness images (flux/mean/median) from the cube using the labels cube.

**shine.shine_utils** — General Utilities
-----------------------------------------

- ``shine_utils.clean_clube``: Continuum subtraction via spectral re-binning and median filtering.
- ``shine_utils.filter_cube``: Spatial and spectral Gaussian smoothing.
- ``shine_utils.subcube``: Spectral sub-cube extraction.
