# -*- coding: utf-8 -*-
__author__ = "Konstantin Klementiev"
__date__ = "16 Feb 2023"

# import os
import sys; sys.path.append('..')  # analysis:ignore
from parseq.core import singletons as csi
from parseq.core import spectra as csp
from . import XAS_nodes as xno
from . import XAS_transforms as xtr
from . import XAS_widgets as xwi


def make_pipeline(withGUI=False):
    csi.withGUI = withGUI

    # instantiate transformation nodes
    nodeIT = xno.NodeIT()
    nodeIF = xno.NodeIF()
    nodeIE = xno.NodeIE()
    nodeMu = xno.NodeMu(xwi.MakeMuWidget if withGUI else None)
    nodeChi = xno.NodeChi(xwi.MakeChiWidget if withGUI else None)
    nodeFT = xno.NodeFT(xwi.MakeFTWidget if withGUI else None)

    # instantiate transformations
    xtr.MakeTrMu(nodeIT, nodeMu)
    xtr.MakeFYMu(nodeIF, nodeMu)
    xtr.MakeTEYMu(nodeIE, nodeMu)
    xtr.MakeChi(nodeMu, nodeChi)
    xtr.MakeFT(nodeChi, nodeFT)

    # initiate data tree
    csi.dataRootItem = csp.Spectrum('root')
    # csi.extraDataFormat['labelName'] = 'label'
