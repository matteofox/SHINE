#!/usr/bin/env python
# coding: utf-8
# AUTHORS: MF, DT
# VERSION: 2.0
#
# SHINE — Spectral Highlighting and Identification of Emission
# =============================================================
# Identifies connected structures in 2D and 3D datasets.
#
# Sub-packages
# -------------
# Find_Em_SHINE : Line-emitter extraction and cataloguing pipeline
# Make_Im_SHINE : 2-D image creation from 3-D cubes

from . import Find_Em_SHINE
from . import Make_Im_SHINE
from .SHINE import runextraction
