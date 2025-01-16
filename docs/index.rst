.. SHINE documentation master file, created by
   sphinx-quickstart on Wed Dec 11 17:33:24 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


==========================
Documentation of SHINE
==========================

**Authors:** Matteo Fossati, Davide Tornotti  

**Date:**  13/12/2024

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
- Dependencies: ``numpy``, ``scipy``, ``astropy``, ``connected-components-3d``

**Steps to Install:**

1. Clone the repository:
   ::

       git clone https://github.com/matteofox/SHINE.git
       cd SHINE

2. Install the code:
   ::

       python -m pip install .


Directory Contents
==================
The github distribution includes a shine/ directory that contains the following codes:

- ``SHINE.py``: The main Python file containing the code for the extraction process.
- ``GUI_SHINE.py``: The Python file containing the code for the GUI.
- ``Make_Im_SHINE.py``: The Python file for the tool used to analyze the extraction results. It generates 2D images.


Usage of SHINE and tools
=================
``SHINE`` and ``Make_Im_SHINE`` are installed as executables that can be called from the terminal.


Run Extraction from the command line
------------------------------------
The extraction is performed using ``SHINE``. The basic idea behind the code is as follows:

- Select a portion of the cube (in the z-direction or wavelength direction) where the user wants to focus the extraction, if needed;
- Mask some voxels based on a user-provided mask (e.g., continuum sources);
- Filter the cube spatially using a user-provided kernel dimension;
- Threshold the cube based on the user-provided S/N threshold;
- Group the connected voxels that satisfy the S/N threshold;
- Write the catalogue and the labeled cube.

**Basic Usage:**

.. code-block:: bash

   SHINE <cube> <variance> --input <input-file> --output <output-file>

**Example (without selecting a subcube):**

.. code-block:: bash

   SHINE ../Cubes/Datacube.fits ../Cubes/Datavar.fits --snthresh 2 --spatsmooth 3 --minvox 300 --mindz 2  --outdir ../Dataproducts/ --writelabels

**Example (selecting a subcube):**

.. code-block:: bash

   SHINE ../Cubes/Datacube.fits ../Cubes/Datavar.fits --zmin 40 --zmax 100 --snthresh 2 --spatsmooth 3 --minvox 300 --mindz 2  --outdir ../Dataproducts/ --writelabels --writesubcube


**Command line arguments for SHINE:**


*General Arguments:*

- ``-h, --help``: Show this help message and exit.

*Input Control Arguments:*

- ``data``: Path of the input data (3D or 2D). Expected to be in extension 0, unless ``extdata`` is defined.
- ``varcube``: Path of the variance cube. Expected to be in extension 0, unless ``extvardata`` is defined.
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

- ``--snthresh``: The SNR of voxels (3D)/pixels (2D) to be included in the extraction (default=2).
- ``--spatsmooth``: Gaussian Sigma of the spatial convolution kernel applied in X and Y (default=0).
- ``--spatsmoothX``: (Optional) Gaussian Sigma of the spatial convolution kernel applied in X. Has priority over ``spatsmooth``.
- ``--spatsmoothY``: (Optional) Gaussian Sigma of the spatial convolution kernel applied in Y. Has priority over ``spatsmooth``.
- ``--specsmooth``: Gaussian Sigma of the spectral convolution kernel applied in Lambda. **(Not implemented yet)**.
- ``--usefftconv``: If ``True``, use FFT for convolution rather than the direct algorithm. For better handling of masked sources, it is currently recommended NOT to use FFT.
- ``--connectivity``: Voxel connectivity scheme to be used (default=26). Allowed values: 4, 8 (2D); 26, 18, 6 (3D).
- ``--maskspedge``: Determines how much (in pixels) to expand the mask around the edges of the cube/image (default=20).

*Cleaning Arguments:*

- ``--minvox``: Minimum number of connected voxels (3D)/pixels (2D) for a source to be in the final catalogue (default=1). For 2D data this argument has priority over minarea.
- ``--mindz``: Minimum number of connected voxels in the spectral direction for a source to be in the final catalogue (default=1). Only valid for 3D data.
- ``--maxdz``: Maximum number of connected voxels in the spectral direction for a source to be in the final catalogue (default=200). Only valid for 3D data.
- ``--minarea``: Minimum number of connected projected spatial voxels (3D)/pixels (2D) for a source to be in the final catalogue (default=1). 

*Output Control Arguments:*

- ``--outdir``: Output directory path.
- ``--writelabels``: If set, write labels cube/image.
- ``--writesmdata``: If set, write the smoothed cube/image.
- ``--writesmvar``: If set, write the smoothed variance.
- ``--writesmsnr``: If set, write the S/N smoothed cube/image.
- ``--writesubcube``: If set and used, write the subcubes (cube and variance). Only valid for 3D data.


Run Extraction using the GUI
----------------------------

Run the GUI using:

.. code-block:: bash

   python GUI_SHINE.py

The GUI is simple (see Fig.1) and allows the user to select input data (Cube/Image, Variance Cube/Image, and optional 2D masks) either by entering the path into the white cells or by clicking the ``Browse`` button to navigate through directories. Similarly, the user can specify the output directory.

All implemented parameters can be adjusted according to the desired extraction, with default values provided for convenience. Data products can be selected by checking the white cells in the *Output Control Arguments* section. By default, the cube containing the IDs of the identified objects and the catalog are selected.

**Note:** Spectral smoothing is not implemented yet.

Once the parameters are set, the user must click on ``Run Script`` to start the extraction. The process takes approximately few minutes for cubes and few seconds for images (faster with FFT convolution; however, for better handling of NaN values, we currently recommend avoiding the use of FFT), after which the output summary is displayed in a new window. The user can then close the GUI and begin analyzing the data products.

.. figure:: /Images/GUI_image.png
   :width: 80%
   :align: left
   :alt: GUI window after ``python GUI_SHINE.py``.

   GUI window after ``python GUI_SHINE.py``.

.. figure:: /Images/Success_GUI_image.png
   :width: 80%
   :align: left
   :alt: Output summary of successful extraction.

   Output summary of successful extraction.

Generation of images
-------------------------------------

This tool is designed to create images from 3D data. It is possible to create three different types of images (<flux>, <mean> or <median>) both selecting only certain voxels based on a 3D mask with associated Ids (e.g. the output labels cube of ``SHINE``, thus creating an extraction image) and all the voxels (thus creating a narrow band image). If <flux> is selected, the units of the output image are :math:`1 \times 10^{-18} \, \mathrm{erg \, s^{-1} \, cm^{-2} \, arcsec^{-2}}`. 
Using a single Id object it is also possible to obtain an image with the representative noise around the object creating a pseudo narrow band around it.
Below is a detailed list of the calling sequence and the command line arguments, grouped by categories:


**Basic Usage:**

.. code-block:: bash

   Make_Im_SHINE <cube> <variance> --input <input-file> --output <output-file>

**Example using a 3D mask for selecting voxels:**

.. code-block:: bash

   Make_Im_SHINE ../Cubes/Datacube.fits ../Cubes/Labels.fits --Id [2,5,9]  --outdir ../Dataproducts/ --itype flux --writeout
   

**Example using a 3D mask for selecting voxels and associate the noise (only for one single object, i.e. one single Id):**

.. code-block:: bash

   Make_Im_SHINE ../Cubes/Datacube.fits ../Cubes/Labels.fits --Id [2]  --outdir ../Dataproducts/ --itype flux --nls -1 --nlsadd 2 --writeout
   

**Example to create a narrow band image:**

.. code-block:: bash

    Make_Im_SHINE ../Cubes/Datacube.fits  --outdir ../Dataproducts/ --itype flux --writeout


**Command line arguments for Make_IM_Shine:**


*Input Control Arguments:*

- ``cube``: Path of the input datacube. Expected to be in extension 0, unless ``extcub`` is defined.
- ``labelsCube``: Path of the cube with labels. Expected to be in extension 0, unless ``extlabels`` is defined.
- ``--Id``: The IDs of the grouped voxels to be used for the surface brightness image extraction. If a list is passed, all the valid IDs will be stacked into the final image. If ``[-1]`` (default), stacks all the IDs in the ``labelsCube``.
- ``--itype``: Type of image to produce. Allowed values: ``flux`` (default), ``mean``, or ``median``.
- ``--extcub``: Specifies the HDU index in the FITS file cube to use for science data extraction (default=0).
- ``--extlabels``: Specifies the HDU index in the FITS file labels to use for labels data extraction (default=0).
- ``--nsl``: Noise layers. Selects this layer associated with the object in the image. If ``-1``, it selects the mean layer. Use this only if ``len(Id)=1``. Default is ``-2`` (no noise layers).
- ``--nsladd``: Noise layers. Specifies how many layers to collapse adjacent to the selected one (default=0).

*Output Control Arguments:*

- ``--outdir``: Output directory path (default=``./``).
- ``--writeout``: If set, writes the flux image and the associated variance image.
- ``--addname``: Optional suffix to append to the base name of the output file. This string will be added after the default filename (e.g., ``flux`` or other predefined parts) and before the file extension. Useful for distinguishing different versions or types of output files.



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

This section describes the available functions in the SHINE module.

.. autofunction:: shine.SHINE.runextraction
   :noindex:

This function performs the extraction process using the specified parameters.

.. autofunction:: shine.SHINE.subcube
   :noindex:

This function extracts a portion of the cube based on user-specified criteria.



License
=======

Copyright (C) 2024 The Authors
  
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.  
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License.

