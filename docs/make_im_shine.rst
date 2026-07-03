Image Generation (``Make_Im_SHINE``)
=====================================

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
