# -*- coding: utf-8 -*-
u"""
Data fits
---------

Linear Combination Fit of µ(E)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: LCF

Function Fit of µ(E)
~~~~~~~~~~~~~~~~~~~~

.. autoclass:: FunctionFit

EXAFS Fit of µ(E)
~~~~~~~~~~~~~~~~~

.. autoclass:: EXAFSFit

"""

__author__ = "Konstantin Klementiev"
__date__ = "15 Feb 2025"
# !!! SEE CODERULES.TXT !!!

import parseq.fits as fits


class LCF(fits.lcf.LCF):
    """
    Work in progress
    """
    name = "LCF"
    ref = "fits.html#linear-combination-fit-of-e"
    tooltip = "Linear Combination Fit. Use a 'flat' normalized mu view."
    xVary = True
    dataAttrs = dict(x='e', y='flat', fit='fitLCF')
    allDataAttrs = dict(x='e', y='flat')
    plotParams = dict(
        fit=dict(linewidth=1.4, linestyle=':', symbol='.', symbolsize=2),
        residue=dict(linewidth=1.0, linestyle='--'))
    nThreads = 4


class FunctionFit(fits.functionfit.FunctionFit):
    """
    Work in progress
    """
    name = "function fit"
    ref = "fits.html#function-fit-of-e"
    tooltip = "Function fit with arbitrary formula. Use a "\
        "'flat' normalized mu view."
    dataAttrs = dict(x='e', y='flat', fit='fitFunc')
    plotParams = dict(fit=dict(linewidth=1., linestyle=':'),
                      residue=dict(linewidth=0.8, linestyle='--'))
    # nThreads = 2
    nProcesses = 2


class EXAFSFit(fits.exafsfit.EXAFSFit):
    """
    Work in progress
    """
    name = "EXAFS fit"
    ref = "fits.html#exafs-fit-of-e"
    tooltip = "EXAFS fit in filtered k- and/or in r-space"
    dataAttrs = dict(x='bftk', y='bft', fit='bftfit',
                     x2='r', y2='ft', fit2='ftfit')
    plotParams = dict(fit=dict(linewidth=1.8, linestyle=':'),
                      # residue=dict(linewidth=1.5, linestyle='--')
                      )
