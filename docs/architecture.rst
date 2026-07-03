Code Architecture
=================

This page provides a visual overview of the SHINE package architecture,
showing how the modules relate to each other and what each one does.


High-Level Architecture
-----------------------

.. graphviz::
   :align: center

   digraph SHINE {
       /* ---- global style ---- */
       graph [fontname="Helvetica", rankdir=TB, nodesep=0.6, ranksep=0.8,
              bgcolor="transparent"];
       node  [fontname="Helvetica", shape=box, style="rounded,filled",
              fontsize=12, margin="0.25,0.12"];
       edge  [fontname="Helvetica", fontsize=10, color="#555555"];

       /* ---- top-level package ---- */
       shine [label="shine\n(package root)", fillcolor="#4A90D9",
              fontcolor=white, fontsize=14, penwidth=2];

       /* ---- core modules ---- */
       core  [label="SHINE.py\nCore extraction engine",
              fillcolor="#5DADE2", fontcolor=white];
       utils [label="shine_utils.py\nUtility functions",
              fillcolor="#5DADE2", fontcolor=white];

       /* ---- sub-packages ---- */
       fem   [label="Find_Em_SHINE\nLine-emitter pipeline",
              fillcolor="#27AE60", fontcolor=white, fontsize=13, penwidth=2];
       mim   [label="Make_Im_SHINE\n2-D image generation",
              fillcolor="#E67E22", fontcolor=white, fontsize=13, penwidth=2];

       /* ---- Find_Em_SHINE modules ---- */
       ext   [label="extraction.py\nStep 1 – SHINE extraction",
              fillcolor="#82E0AA"];
       cov   [label="covariance.py\nStep 2 – Noise covariance",
              fillcolor="#82E0AA"];
       cat   [label="catalogue.py\nStep 3 – Catalogue & cutouts",
              fillcolor="#82E0AA"];
       cli   [label="cli.py\nCommand-line interface",
              fillcolor="#82E0AA"];

       /* ---- Make_Im_SHINE modules ---- */
       mmod  [label="Make_Im_SHINE.py\nImage creation",
              fillcolor="#F0B27A"];

       /* ---- edges ---- */
       shine -> core  [label=" core"];
       shine -> utils [label=" utils"];
       shine -> fem   [label=" sub-pkg"];
       shine -> mim   [label=" sub-pkg"];

       fem -> ext [label=" step 1"];
       fem -> cov [label=" step 2"];
       fem -> cat [label=" step 3"];
       fem -> cli [label=" CLI"];

       mim -> mmod;

       /* ---- cross-module dependencies ---- */
       edge [style=dashed, color="#999999"];
       ext  -> core  [label="calls runextraction"];
       ext  -> utils [label="calls filter_cube,\nclean_clube"];
       cat  -> cov   [label="uses covariance\ncorrection"];
   }


The diagram above shows the package tree with the main call dependencies.
Solid arrows indicate ownership (parent → child module); dashed arrows
indicate cross-module function calls.


Module Summary
--------------

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Component
     - Module
     - Purpose
   * - **Core**
     - ``SHINE.py``
     - Extraction engine: filtering, thresholding, connected-component
       labelling, catalogue generation, photometry. Entry point:
       :func:`~shine.SHINE.runextraction`.
   * - **Utilities**
     - ``shine_utils.py``
     - General-purpose helpers: continuum subtraction
       (``clean_clube``), spatial/spectral smoothing (``filter_cube``),
       sub-cube cropping (``subcube``).
   * - **Find_Em_SHINE**
     - ``extraction.py``
     - *Step 1* — Wrapper around ``runextraction`` with emitter-optimised
       defaults; optional continuum subtraction.
   * -
     - ``covariance.py``
     - *Step 2* — Monte Carlo estimation of the empirical noise covariance
       as a function of aperture size and wavelength.
   * -
     - ``catalogue.py``
     - *Step 3* — S/N correction, quality cuts, confidence-class assignment,
       per-source image cutouts and spectral extraction.
   * -
     - ``cli.py``
     - Command-line interface that orchestrates the three steps from a
       single ``config.ini`` file.
   * - **Make_Im_SHINE**
     - ``Make_Im_SHINE.py``
     - Generates 2-D surface-brightness images from cubes and
       segmentation maps.


Find_Em_SHINE Pipeline Flow
----------------------------

The emitter-detection pipeline is designed to run sequentially in three
self-contained steps.  Each step produces output files that are consumed
by the next:

.. graphviz::
   :align: center

   digraph pipeline {
       graph [fontname="Helvetica", rankdir=LR, nodesep=0.5, ranksep=1.0,
              bgcolor="transparent"];
       node  [fontname="Helvetica", shape=box, style="rounded,filled",
              fontsize=11, margin="0.2,0.1"];
       edge  [fontname="Helvetica", fontsize=9];

       /* ---- optional pre-processing ---- */
       pre  [label="Pre-processing\n(optional)\nclean_clube\nfilter_cube",
             fillcolor="#D5F5E3", style="rounded,filled,dashed"];

       /* ---- pipeline steps ---- */
       s1   [label="Step 1\nextraction.py\n──────────\nSegmentation map\nFiltered cubes\nRaw catalogue",
             fillcolor="#82E0AA"];
       s2   [label="Step 2\ncovariance.py\n──────────\nCovariance .npz\nPolynomial fit\nDiagnostic plot",
             fillcolor="#82E0AA"];
       s3   [label="Step 3\ncatalogue.py\n──────────\nS/N-corrected cat.\nConfidence classes\nImage cutouts\n1-D spectra",
             fillcolor="#82E0AA"];

       pre -> s1 [label="cubes"];
       s1  -> s2 [label="filtered cubes\n+ labels"];
       s1  -> s3 [label="raw catalogue\n+ seg. map"];
       s2  -> s3 [label="covariance\npolynomial"];
   }


File Tree
---------

For reference, the on-disk layout of the ``shine/`` package::

    shine/
    ├── __init__.py               # Package root
    ├── SHINE.py                  # Core extraction engine
    ├── shine_utils.py            # Utility functions
    │
    ├── Find_Em_SHINE/            # Line-emitter pipeline
    │   ├── __init__.py
    │   ├── cli.py                # CLI entry-point
    │   ├── extraction.py         # Step 1 – extraction
    │   ├── covariance.py         # Step 2 – covariance
    │   └── catalogue.py          # Step 3 – catalogue
    │
    └── Make_Im_SHINE/            # 2-D image generation
        ├── __init__.py
        └── Make_Im_SHINE.py
