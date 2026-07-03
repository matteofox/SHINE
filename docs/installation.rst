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
