#!/usr/bin/env python
# coding: utf-8
# AUTHORS: MF, DT
# VERSION: 2.0
#
# shine_utils.py
# --------------
# General-purpose utility functions for the SHINE pipeline.
#
# Functions specific to the emitter-extraction workflow (covariance
# estimation, catalogue building, etc.) have been moved to the
# ``shine.Find_Em_SHINE`` sub-package.


import warnings

import numpy as np
from scipy.ndimage import median_filter

from astropy.stats import sigma_clipped_stats

warnings.filterwarnings("ignore", category=RuntimeWarning)


# =============================================================================
# 1.  clean_clube — continuum subtraction
# =============================================================================

def clean_clube(data, filtsize=7, rebinfac=40):
    """Subtract a smooth continuum from a 3-D spectroscopic cube.

    The cube is first re-binned spectrally by collapsing groups of
    ``rebinfac`` layers into a single median-stack slice. The resulting
    low-resolution continuum cube is then smoothed along the spectral axis
    with a sliding-window median filter of width ``filtsize``. The smoothed
    continuum is finally subtracted layer-by-layer from the original cube.

    Parameters
    ----------
    data : numpy.ndarray, shape (nz, ny, nx)
        Input 3-D data cube.  NaN values are treated as masked pixels via
        :class:`numpy.ma.MaskedArray`.
    filtsize : int, optional
        Width (in re-binned spectral slices) of the median filter applied to
        the low-resolution continuum cube. Default is 7.
    rebinfac : int, optional
        Number of original spectral layers collapsed into each continuum
        slice. Default is 40.

    Returns
    -------
    data : numpy.ma.MaskedArray, shape (nz, ny, nx)
        Continuum-subtracted cube.  The mask mirrors the NaN pattern of the
        input.
    """
    nz, ny, nx = np.shape(data)

    data = np.ma.array(data, mask=np.isnan(data))

    zrebin = int(np.ceil(nz / rebinfac))

    contcube = np.zeros((zrebin, ny, nx))

    print(f'... Rebinning the cube into {zrebin} continuum slices '
          f'(rebinfac={rebinfac})')
    for ii in np.arange(zrebin):
        zmin = rebinfac * ii
        zmax = min(nz, rebinfac * (ii + 1))
        _, median, _ = sigma_clipped_stats(
            data[zmin:zmax, :, :], sigma=3, axis=0, maxiters=3
        )
        contcube[ii] = median

    print(f'... Filtering the continuum cube with a spectral median filter '
          f'(filtsize={filtsize})')
    filtcube = median_filter(contcube, size=filtsize, axes=0)

    print('... Subtracting the smooth continuum')
    for ii in np.arange(zrebin):
        zmin = rebinfac * ii
        zmax = min(nz, rebinfac * (ii + 1))
        data[zmin:zmax, :, :] -= filtcube[ii]

    return data
