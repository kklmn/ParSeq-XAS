# -*- coding: utf-8 -*-
"""
A pipeline for data processing of XAS spectra.
"""

import os.path as osp

import sys; sys.path.append('..')  # analysis:ignore
from parseq.core import singletons as csi
from .XAS_pipeline import make_pipeline
from .XAS_tests import load_test_data

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
