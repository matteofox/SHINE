.. SHINE documentation master file, created by
   sphinx-quickstart on Wed Dec 11 17:33:24 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


==========================
Documentation of SHINE
==========================

**Authors:** Matteo Fossati, Davide Tornotti  
**Date:**  

.. contents::
   :depth: 2
   :local:


Introduction
============

**Project Name:** ``SHINE``  
SHINE (Spectral Highlighting and Identification of Emission) is a simple Python-based script that identifies connected structures above a user-given signal-to-noise (S/N) threshold in 2D and 3D datasets.


Installation
============

**Requirements:**

- Python version: ``>=3.8``
- Dependencies: ``numpy``, ``scipy``, ``astropy``, ``cc3d``

**Steps to Install:**

1. Clone the repository:
   ::

       git clone https://github.com/matteofox/SHINE.git
       cd SHINE

2. Install dependencies:
   ::

       pip install -r requirements.txt


Directory Contents
==================

- ``SHINE.py``: The main Python file containing the code for the extraction process.
- ``GUI_SHINE.py``: The Python file containing the code for the GUI.
- ``Make_Im_SHINE.py``: The Python file for the initial tool used to analyze the extraction results. It generates 2D surface brightness images.


Usage of SHINE.py
=================

**Basic Usage:**

.. code-block:: bash

   python SHINE.py <cube> <variance> --input <input-file> --output <output-file>

**Example:**

.. code-block:: bash

   python SHINE.py ../Cubes/Datacube.fits ../Cubes/Datavar.fits --snthresh 2 --spatsmooth 3 --minvox 300 --mindz 2  --outdir ../Dataproducts/ --writelabels

**GUI Usage:**

Run the GUI using:

.. code-block:: bash

   python GUI_SHINE.py

The GUI is simple and allows the user to select input data (Cube, Variance, and optional 2D masks) either by entering the path into the white cells or by clicking the ``Browse`` button to navigate through directories. Similarly, the user can specify the output directory.

All implemented parameters can be adjusted according to the desired extraction, with default values provided for convenience. Data products can be selected by checking the white cells in the *Output Control Arguments* section. By default, the cube containing the IDs of the identified objects and the catalog are selected.

**Note:** Spectral smoothing is not implemented yet.

Once the parameters are set, the user must click on ``Run Script`` to start the extraction. The process takes approximately one minute (faster with FFT convolution; however, for better handling of NaN values, we currently recommend avoiding the use of FFT), after which the output summary is displayed in a new window. The user can then close the GUI and begin analyzing the data products.


Features of SHINE.py
=====================

This tool provides several features and arguments for controlling its behavior. Below is a detailed list of the available arguments, grouped by their respective categories:

**General Arguments:**

- ``-h, --help``: Show this help message and exit.

**Input Control Arguments:**

- ``cube``: Path of the input datacube. Expected to be in extension 0, unless ``extcub`` is defined.
- ``varcube``: Path of the variance cube. Expected to be in extension 0, unless ``extvar`` is defined.
- ``--mask2d``: Path of an optional two-dimensional mask to be applied along the wave axis.
- ``--mask2dpost``: Path of an optional two-dimensional mask to be applied after smoothing along the wave axis.
- ``--mask3d``: Path of an optional three-dimensional mask. **(Not implemented yet)**.
- ``--extcub``: Specifies the HDU index in the FITS file cube to use for data cube extraction.
- ``--extvar``: Specifies the HDU index in the FITS file variance to use for cube extraction.

**Extraction Arguments:**

- ``--snthresh``: The SNR of voxels to be included in the extraction.
- ``--spatsmooth``: Gaussian Sigma of the spatial convolution kernel applied in X and Y.
- ``--spatsmoothX``: Gaussian Sigma of the spatial convolution kernel applied in X. Has priority over ``spatsmooth``.
- ``--spatsmoothY``: Gaussian Sigma of the spatial convolution kernel applied in Y. Has priority over ``spatsmooth``.
- ``--specsmooth``: Gaussian Sigma of the spectral convolution kernel applied in Lambda. **(Not implemented yet)**.
- ``--usefftconv``: If ``True``, use FFT for convolution rather than the direct algorithm.
- ``--connectivity``: Voxel connectivity scheme to be used. Allowed values: 4, 8 (2D); 26, 18, 6 (3D).
- ``--maskspedge``: Determines how much (in pixels) to expand the mask around the edges of the cube/image.

**Cleaning Arguments:**

- ``--minvox``: Minimum number of connected voxels for a source to be in the final catalogue.
- ``--mindz``: Minimum number of connected voxels in the spectral direction for a source to be in the final catalogue.
- ``--maxdz``: Maximum number of connected voxels in the spectral direction for a source to be in the final catalogue.
- ``--minarea``: Minimum number of connected voxels in the projected spatial direction for a source to be in the final catalogue.

**Output Control Arguments:**

- ``--outdir``: Output directory path.
- ``--writelabels``: If set, write labels cube.
- ``--writesmcube``: If set, write the smoothed cube.
- ``--writesmvar``: If set, write the smoothed variance.
- ``--writesmsnrcube``: If set, write the S/N smoothed cube.


Usage of Make_Im_SHINE.py
=========================

**Basic Usage:**

.. code-block:: bash

   python Make_Im_SHINE.py <cube> <variance> --input <input-file> --output <output-file>

**Example:**

.. code-block:: bash

   python Make_Im_SHINE.py ../Cubes/Datacube.fits ../Cubes/Datavar.fits --Id [2,5,9]  --outdir ../Dataproducts/ --writeout


Features of Make_Im_SHINE.py
============================

This tool is designed to create surface brightness images from 3D data using the output of ``SHINE``. The units of the output image are :math:`1 \times 10^{-18} \, \mathrm{erg \, s^{-1} \, cm^{-2} \, arcsec^{-2}}`. Below is a detailed list of the available arguments, grouped by their respective categories:

**Input Control Arguments:**

- ``cube``: Path of the input datacube. Expected to be in extension 0, unless ``extcub`` is defined.
- ``varcube``: Path of the variance cube. Expected to be in extension 0, unless ``extvar`` is defined.
- ``labelsCube``: Path of the cube containing the labels. Expected to be in extension 0, unless ``extlabels`` is defined.
- ``--Id``: IDs of the grouped voxels to be used for surface brightness image extraction. If a list of IDs is provided, all valid IDs will be stacked into the final image. Default is ``[-1]``.
- ``--extcub``: Specifies the HDU index in the FITS file for the science data. Default is ``0``.
- ``--extvar``: Specifies the HDU index in the FITS file for the variance data. Default is ``0``.
- ``--extlabels``: Specifies the HDU index in the FITS file for the labels data. Default is ``0``.

**Output Control Arguments:**

- ``--outdir``: Specifies the output directory path. Default is ``./``.
- ``--writeout``: If set, the tool will write the flux image and the associated variance image.


Contributing
============

Explain how others can contribute to your project:

1. Fork the repository on GitHub.
2. Create a new branch for your feature/bugfix.
3. Submit a pull request.


License
=======

Copyright (C) 2024 The Authors  
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.  
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public

