# -*- coding: utf-8 -*-
__author__ = "Konstantin Klementiev"
__date__ = "20 Nov 2023"

import sys; sys.path.append('..')  # analysis:ignore
from parseq.core import singletons as csi
from parseq.core import spectra as csp

from parseq.gui.fits.glcf import LCFWidget
from parseq.gui.fits.gfunctionfit import FunctionFitWidget
from parseq.gui.fits.gexafsfit import EXAFSFitWidget

from . import XAS_nodes as xno
from . import XAS_transforms as xtr
from . import XAS_widgets as xwi
from . import XAS_fits as xfi


def make_pipeline(withGUI=False):
    csi.withGUI = withGUI

    # instantiate transformation nodes
    nodeIT = xno.NodeIT(xwi.CurWidget if withGUI else None)
    nodeIF = xno.NodeIF(xwi.CurWidget if withGUI else None)
    nodeIE = xno.NodeIE(xwi.CurWidget if withGUI else None)
    nodeIXES = xno.NodeIXES(xwi.HERFDWidget if withGUI else None)
    nodeMu = xno.NodeMu(
        [xwi.MuWidget, xwi.MuSelfAbsorptionCorrection] if withGUI else None)
    nodeChi = xno.NodeChi(xwi.ChiWidget if withGUI else None)
    nodeFT = xno.NodeFT(xwi.FTWidget if withGUI else None)
    nodeBFT = xno.NodeBFT()

    # instantiate transformations
    xtr.MakeTrMu(nodeIT, nodeMu)
    xtr.MakeFYMu(nodeIF, nodeMu)
    xtr.MakeTEYMu(nodeIE, nodeMu)
    xtr.MakeHERFD(nodeIXES, nodeMu)
    xtr.MakeChi(nodeMu, nodeChi)
    xtr.MakeFT(nodeChi, nodeFT)
    xtr.MakeBFT(nodeFT, nodeBFT)

    # instantiate fits
    xfi.LCF_mu(nodeMu, LCFWidget if withGUI else None)
    xfi.FunctionFit_mu(nodeMu, FunctionFitWidget if withGUI else None)
    xfi.EXAFSFit([nodeBFT, nodeFT], EXAFSFitWidget if withGUI else None)

    # initiate data tree
    csi.dataRootItem = csp.Spectrum('root')
