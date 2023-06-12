# -*- coding: utf-8 -*-
__author__ = "Konstantin Klementiev"
__date__ = "24 May 2023"
# !!! SEE CODERULES.TXT !!!

import parseq.fits as fits


class LCF(fits.lcf.LCF):
    name = 'LCF'
    tooltip = "Linear Combination Fit. Use a 'flat' normalized mu view."
    xVary = True
    dataAttrs = dict(x='e', y='flat', fit='fitLCF')
    allDataAttrs = dict(x='e', y='flat')
    plotParams = dict(
        fit=dict(linewidth=1.4, linestyle=':', symbol='.', symbolsize=2),
        residue=dict(linewidth=1.0, linestyle='--'))
    nThreads = 4


class FunctionFit(fits.functionfit.FunctionFit):
    name = 'function fit'
    tooltip = "Function fit with arbitrary formula. Use a "\
        "'flat' normalized mu view."
    dataAttrs = dict(x='e', y='flat', fit='fitFunc')
    plotParams = dict(fit=dict(linewidth=1., linestyle=':'),
                      residue=dict(linewidth=0.8, linestyle='--'))
    nThreads = 2