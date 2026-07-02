#!/usr/bin/env python
# coding: utf-8
# AUTHORS: MF, DT
# VERSION: 2.0
#
# Find_Em_SHINE sub-package
# -------------------------
# End-to-end tools for the extraction and cataloguing of line emitters
# from spectroscopic cubes using the SHINE pipeline.
#
# Three main steps:
#   1. extract()          — run SHINE extraction with emitter-optimised defaults
#   2. covariance()       — estimate the empirical noise covariance
#   3. build_em_catalog() — build the final emitter catalogue with corrected S/N

from .extraction import extract
from .covariance import estimate_empirical_covariance as covariance
from .catalogue import build_emitter_catalogue as build_em_catalog
