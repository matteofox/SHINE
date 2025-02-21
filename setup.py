from setuptools import setup, find_packages

setup(
    name='SHINE',  # Replace with your project name
    version='1.1',             # Initial version
    packages=find_packages(),  # Automatically finds all the sub-packages in the project
    install_requires=[         # List of dependencies
        'numpy',  # Example dependency
        'astropy',
        'scipy',
        'connected-components-3d',  # Example dependency
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU GPL 3',  # Specify your license
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',  # Minimum Python version required
    
    entry_points={
        'console_scripts': [
            'SHINE = shine.SHINE:main',  # If you have CLI commands
            'Make_Im_SHINE = shine.Make_Im_SHINE:main'
        ],
    },
    author='Matteo Fossati, Davide Tornotti',  # Replace with your name
    author_email='matteo.fossati@unimib.it',  # Replace with your email
    description='Spectral Highlighting and Identification of Emission identifies connected structures in 2D and 3D datasets.',
    long_description=open('README.rst').read(),  # Optional: use a README file for long description
    long_description_content_type='text/markdown',  # Assuming README.md is in markdown format
    url='https://github.com/matteofox/SHINE',  # Replace with your project URL
)
