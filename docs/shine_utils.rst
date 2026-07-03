General Utilities (``shine_utils``)
=====================================

The ``shine_utils`` module provides general-purpose utility functions for spectroscopic cubes.


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
