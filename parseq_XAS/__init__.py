# -*- coding: utf-8 -*-
"""
The pipeline calculates absorption coefficient from transmission signals or
signals of secondary processes (fluorescence or electron yield), determines
absorption edge position, calculates pre- and post-edge background, calculates
EXAFS function, its forward and backward Fourier Transform. The pipeline has
a Linear Combination Fit and an EXAFS fit.
"""

import os.path as osp

import sys; sys.path.append('..')  # analysis:ignore
from parseq.core import singletons as csi

from .version import __versioninfo__, __version__, __date__

__author__ = "Konstantin Klementiev (MAX IV Laboratory)"
__email__ = "first dot last at gmail dot com"
__license__ = "MIT license"
__synopsis__ = "A pipeline for data processing of XAS spectra"

csi.pipelineName = 'XAS'
csi.appPath = osp.dirname(osp.abspath(__file__))
csi.appIconPath = osp.join(
    csi.appPath, 'doc', '_images', 'XAS_icon.ico')
csi.appSynopsis = __synopsis__
csi.appDescription = __doc__
csi.appAuthor = __author__
csi.appLicense = __license__
csi.appVersion = __version__

from .XAS_pipeline import make_pipeline
from .XAS_tests import load_test_data
