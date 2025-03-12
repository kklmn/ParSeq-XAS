# -*- coding: utf-8 -*-
u"""
No GUI: data fits
-----------------

Linear Combination Fit of µ(E)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: LCF_mu

Function Fit of µ(E)
~~~~~~~~~~~~~~~~~~~~

.. autoclass:: FunctionFit_mu

EXAFS Fit of χ(k) and χ(r)
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: EXAFSFit

"""

__author__ = "Konstantin Klementiev"
__date__ = "15 Feb 2025"
# !!! SEE CODERULES.TXT !!!

import parseq.fits as fits


class LCF_mu(fits.lcf.LCF):
    """
    The Linear Combination Fit (LCF) fits a weighted sum of reference spectra
    to the spectrum of interest. Here, the array attribute :param:`flat` (the
    version of µ that has a horizontal post-edge) of the involved data objects
    are used as ordinates and the attribute :param:`e` is used as abscissa. The
    fit curve resides in a new array :param:`fitLCF` as an attribute of the
    data container.

    The fit creates an attribute :param:`lcf_params` of the data container that
    is a list of dictionaries, a dictionary per reference spectrum, with the
    following items:
       | *name*: str, the alias of the reference spectrum,
       | *use*: bool, switches the reference on/off,
       | *w*: float, the weight of the reference,
       | *wBounds*: 3-list, [*min*, *max*, *Δ*], Δ is only used by the GUI,
       |
       | If :param:`xVary` is True (see below), these two entries are needed:
       | *dx*: float, the energy shift of the reference,
       | *dxBounds*: 3-list, [*min*, *max*, *Δ*] for the energy shift.

    The reference dictionary may have *wtie* and *dxtie* for tie expressions
    for the weight *w* and the energy shift *dx*. Tie expression is a string
    variable either starting with "fix" or starting with one of the three signs
    "=><" and followed by a Python expression of other weights accessed as
    "w[1]", "w[2]", "dx[1]" etc. It is possible to add metavariables with a
    user-defined name to be used in tie expressions. A variable can also be
    fixed by its *min* and *max* bounds: if *min* >= *max*, the value will be
    fixed at the value of *min*.

    After the fit, the items *wError* and *dxError* contain the fitting errors.

    The class variable :param:`xVary` determines whether the reference spectra
    are allowed to vary their energy axis. This feature may be needed if the
    energy calibration of the involved spectra is not trustful.

    The actual curve fitting is done by ``scipy.optimize.curve_fit()``. The
    attribute :param:`lcf_result` of the data container will get the fit info
    returned by ``curve_fit()``.
    """

    name = "LCF"
    ref = "nogui.html#linear-combination-fit-of-e"
    tooltip = "Linear Combination Fit. Use a 'flat' normalized mu view."
    xVary = True
    dataAttrs = dict(x='e', y='flat', fit='fitLCF')
    allDataAttrs = dict(x='e', y='flat')
    plotParams = dict(
        fit=dict(linewidth=1.4, linestyle=':', symbol='.', symbolsize=2),
        residue=dict(linewidth=1.0, linestyle='--'))
    nThreads = 4


class FunctionFit_mu(fits.functionfit.FunctionFit):
    """
    The Function Fit fits a Python expression that depends on a set of user
    variables to data. Here, the array attribute :param:`flat` (the version of
    µ that has a horizontal post-edge) of the fitted data object is used as
    ordinate and the attribute :param:`e` is used as abscissa. The fit curve
    resides in a new array :param:`fitFunc` as an attribute of the data
    container.

    The fit creates a few attributes of the data container:

        * :param:`ffit_formula`: str, fit formula,
        * :param:`ffit_params`: dictionary of fitting variables as keys,
        * :param:`ffit_xRange`: 2-list [*min*, *max*] of energy limits.

    The dictionary :param:`ffit_params` is a dictionary of dictionaries of the
    following keys: *value* (float, the fit variable value), *lim* (2-list of
    *min* and *max*), *tie* (str, a tie expression) and *error* (float, the
    fitting error).

    The actual curve fitting is done by ``scipy.optimize.curve_fit()``. The
    attribute :param:`ffit_result` of the data container will get the fit info
    returned by ``curve_fit()``.
    """

    name = "function fit"
    ref = "nogui.html#function-fit-of-e"
    tooltip = "Function fit with arbitrary formula. Use a "\
        "'flat' normalized mu view."
    dataAttrs = dict(x='e', y='flat', fit='fitFunc')
    plotParams = dict(fit=dict(linewidth=1., linestyle=':'),
                      residue=dict(linewidth=0.8, linestyle='--'))
    # nThreads = 2
    nProcesses = 2


class EXAFSFit(fits.exafsfit.EXAFSFit):
    """
    The EXAFS Fit fits a sum of EXAFS shells to data. The array attribute
    :param:`bft` of the fitted data object is used as ordinate and the
    attribute :param:`bftk` is used as abscissa. The fit curve resides in a new
    array :param:`bftfit` as an attribute of the data container. Optionally,
    the fit can be done in the Fourier space; then the fit works with
    :param:`ft` vs :param:`r` and the fit curve is :param:`ftfit`.

    The fit creates a few attributes of the data container:

        * :param:`exafsfit_params`: list of shell dictionaries,
        * :param:`exafsfit_aux`: list of shell specifications -- a path to
          the corresponding FEFF file and its descriptions,
        * :param:`exafsfit_kRange`: 2-list [*min*, *max*]
        * :param:`exafsfit_k_use`: bool, whether to fit in k-space
        * :param:`exafsfit_rRange`: 2-list [*min*, *max*]
        * :param:`exafsfit_r_use`: bool, whether to fit in r-space

    The list :param:`exafsfit_params` contains shell dictionaries, per shell.
    Each shell dictionary defines the 4 shell fitting parameters 'r', 'n', 's',
    'e' as keys and fitting parameter dictionaries as values. The fitting
    parameter dictionary defines the pairs 'value': float, 'step': float,
    'lim': [min, max] and 'error': float. It may additionally have a
    ('tie': expression) entry.

    The list :param:`exafsfit_params` has an extra shell dictionary that
    contains an 's0' fitting variable (multi-electron damping factor) and
    optional metavariables.

    The actual curve fitting is done by ``scipy.optimize.curve_fit()``. The
    attribute :param:`exafsfit_result` of the data container will get the fit
    info returned by ``curve_fit()``.
    """

    name = "EXAFS fit"
    ref = "nogui.html#exafs-fit-of-k-and-r"
    tooltip = "EXAFS fit in filtered k- and/or in r-space"
    dataAttrs = dict(x='bftk', y='bft', fit='bftfit',
                     x2='r', y2='ft', fit2='ftfit')
    plotParams = dict(fit=dict(linewidth=1.8, linestyle=':'),
                      # residue=dict(linewidth=1.5, linestyle='--')
                      )
