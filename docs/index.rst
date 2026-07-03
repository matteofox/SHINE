.. SHINE documentation master file.
   Authors: Matteo Fossati, Davide Tornotti

==========================
Documentation of SHINE
==========================

**Authors:** Matteo Fossati, Davide Tornotti

**Date:** 03/07/2026 (v2.1)

.. warning::
   **Version in Development (Unstable Branch)**

   This documentation refers to a version of the code currently under development, restructuring and testing (non-main branch). Some features may not be stable or verified yet.

   For the stable, tested and production-ready version, please refer to the ``main`` branch and the official documentation at `this link <https://shinespec.readthedocs.io/en/latest/>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   shine_core
   find_em_shine
   make_im_shine
   shine_utils
   api
   changelog


Introduction
============

**Project Name:** ``SHINE``

``SHINE`` (Spectral Highlighting and Identification of Emission) is a Python package that identifies connected structures above a user-given signal-to-noise (S/N) threshold in 2D and 3D datasets, and provides a full suite of data-analysis utilities for spectroscopic cubes.

Version 2.1 introduces a modular sub-package architecture:

* **Core SHINE** — Core extraction engine for 2D and 3D data (filtering,
  thresholding, connected-component labelling, catalogue generation).
* **Find_Em_SHINE** — End-to-end pipeline for line-emitter detection and
  cataloguing, organised in three steps: extraction, covariance estimation,
  and catalogue building.
* **Make_Im_SHINE** — Tool to generate 2-D surface-brightness images from
  spectroscopic cubes using the SHINE extraction products.

Previous changes (v2.0):

* **shine_utils** — general-purpose utility functions (continuum subtraction, spatial/spectral smoothing, sub-cube cropping).
* Removal of the Tkinter GUI (``GUI_SHINE``).
* ``matplotlib`` added as a core dependency; ``sep`` added as an optional
  dependency for automatic source masking in the covariance estimator.
* Fully revised NumPy-style docstrings for all public functions.


Package Structure
==================

The ``shine/`` package is organised as follows::

    shine/
    ├── __init__.py             # Package root: exposes sub-packages and runextraction
    ├── SHINE.py                # Core extraction engine (runextraction)
    ├── shine_utils.py          # General-purpose utilities (clean_clube, filter_cube, subcube)
    │
    ├── Find_Em_SHINE/          # Line-emitter extraction and cataloguing
    │   ├── __init__.py         # Exposes: extract, covariance, build_em_catalog
    │   ├── cli.py              # CLI for emitters with config.ini
    │   ├── extraction.py       # Step 1: SHINE extraction with emitter defaults
    │   ├── covariance.py       # Step 2: empirical noise-covariance estimation
    │   └── catalogue.py        # Step 3: final catalogue, S/N correction, cutouts
    │
    └── Make_Im_SHINE/          # 2-D image generation
        ├── __init__.py         # Exposes: Make_Im_SHINE, main
        └── Make_Im_SHINE.py    # Image creation from cubes and labels


Contributing
============

If you are interested in contributing to the project, please contact us and follow these steps:

1. Fork the repository on GitHub.
2. Create a new branch for your feature/bugfix.
3. Submit a pull request.


License
=======

Copyright (C) 2024–2026 The Authors

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License.
