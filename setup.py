from setuptools import setup, find_packages

setup(
    name='SHINE',
    version='2.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'astropy',
        'scipy',
        'matplotlib',
        'connected-components-3d',
    ],
    extras_require={
        # sep is used by shine_utils.estimate_empirical_covariance for
        # automatic continuum-source detection.  Install it with:
        #   pip install sep
        # If sep is not installed, you must provide a mask_source explicitly.
        'covariance': ['sep'],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',
    entry_points={
        'console_scripts': [
            'SHINE = shine.SHINE:main',
            'Make_Im_SHINE = shine.Make_Im_SHINE.Make_Im_SHINE:main',
        ],
    },

    author='Matteo Fossati, Davide Tornotti',
    author_email='matteo.fossati@unimib.it',
    description=(
        'Spectral Highlighting and Identification of Emission — identifies '
        'connected structures in 2D and 3D datasets.'
    ),
    long_description=open('README.rst').read(),
    long_description_content_type='text/x-rst',
    url='https://github.com/matteofox/SHINE',
)

