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
__date__ = "13 Mar 2025"
# !!! SEE CODERULES.TXT !!!

import parseq.fits as fits


class LCF_mu(fits.lcf.LCF):
    """
    The Linear Combination Fit (LCF) models the spectrum of interest as a
    weighted sum of reference spectra. In this procedure, the array attribute
    :param:`flat` (the flattened absorption spectrum with a horizontal
    post-edge) is used as the ordinate, while :param:`e` serves as the
    abscissa. The resulting fit curve is stored in a new array attribute,
    :param:`fitLCF`, within the data container.

    The fitting process creates an attribute :param:`lcf_params` in the data
    container. This attribute is a list of dictionaries -- one per reference
    spectrum -- with the following entries:

    *name* str, alias of the reference spectrum

    *use* bool, enables or disables the reference

    *w* float, weight of the reference

    *wBounds* 3-list [*min*, *max*, *Δ*]; Δ is used only by the GUI

    If :param:`xVary` is enabled, the following additional entries are
    included:

    *dx* float, energy shift of the reference

    *dxBounds* 3-list [*min*, *max*, *Δ*] for the energy shift

    Optional entries *wtie* and *dxtie* define tie expressions for the weight
    (*w*) and energy shift (*dx*). A tie expression is a string that either
    begins with "fix" or with one of the symbols "=", ">", or "<", followed by
    a Python expression involving other parameters (e.g. "w[1]", "dx[2]").
    User-defined metavariables may also be introduced for use in such
    expressions. Alternatively, parameters can be fixed by setting their bounds
    such that *min* ≥ *max*, in which case the parameter is fixed at the value
    of *min*.

    After fitting, the entries *wError* and *dxError* store the estimated
    uncertainties of the corresponding parameters.

    The class variable :param:`xVary` determines whether reference spectra are
    allowed to shift along the energy axis, which can be useful when energy
    calibration is uncertain.

    Optionally, a pinhole fraction can be introduced as an additional fitting
    parameter, primarily relevant for transmission measurements.

    The fitting itself is performed using ``scipy.optimize.curve_fit()``. The
    resulting fit information returned by this function is stored in the
    :attr:`lcf_result` attribute of the data container.
    """

    name = "LCF"
    ref = "nogui.html#linear-combination-fit-of-e"
    tooltip = "Linear Combination Fit. Use a 'flat' normalized mu view."
    xVary = True
    dataAttrs = dict(x='e', y='flat', fit='fitLCF',
                     pre_edge='pre_edge', post_edge='post_edge')
    allDataAttrs = dict(x='e', y='flat')
    plotParams = dict(
        fit=dict(linewidth=1.4, linestyle=':', symbol='.', symbolsize=2),
        residue=dict(linewidth=1.0, linestyle='--'))
    nThreads = 4


class FunctionFit_mu(fits.functionfit.FunctionFit):
    """
    The Function Fit fits a Python expression that depends on a set of user
    variables to the data. In this procedure, the array attribute
    :param:`flat` (the flattened absorption coefficient with a horizontal
    post-edge) of the fitted data object is used as the ordinate, while
    :param:`e` serves as the abscissa. The fit curve is stored as a new array
    attribute, :param:`fitFunc`, in the data container.

    The fitting process creates several attributes in the data container:

    * :param:`ffit_formula` str, the fit formula

    * :param:`ffit_params` dict, fitting variables as keys

    * :param:`ffit_xRange` 2-list [*min*, *max*], defining the energy range

    The dictionary :param:`ffit_params` is itself a dictionary of dictionaries
    with the following entries:

    *value* float, current value of the fitting variable

    *lim* 2-list [*min*, *max*], bounds of the variable

    *tie* str, tie expression linking the variable to others

    *error* float, estimated uncertainty of the variable

    The fitting itself is performed using ``scipy.optimize.curve_fit()``. The
    resulting fit information returned by this function is stored in the
    :param:`ffit_result` attribute of the data container.
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
    The EXAFS Fit models the data as a sum of EXAFS shells. The array
    attribute :param:`bft` of the fitted data object is used as the ordinate,
    while the attribute :param:`bftk` serves as the abscissa. The fit curve is
    stored as a new array attribute, :param:`bftfit`, in the data container.
    Optionally, the fit can be performed in Fourier space; in this case,
    :param:`ft` is used versus :param:`r`, and the fit curve is stored in
    :param:`ftfit`.

    The fitting process creates several attributes in the data container:

    * :param:`exafsfit_params` list of dictionaries, one per EXAFS shell

    * :param:`exafsfit_aux` list of shell specifications, including paths to
      FEFF files and their descriptions

    * :param:`exafsfit_kRange` 2-list [*min*, *max*] of the k-range

    * :param:`exafsfit_k_use` bool, whether the fit is performed in k-space

    * :param:`exafsfit_rRange` 2-list [*min*, *max*] of the r-range

    * :param:`exafsfit_r_use` bool, whether the fit is performed in r-space

    The list :param:`exafsfit_params` contains one dictionary per shell. Each
    shell dictionary defines the four fitting parameters ``r``, ``n``, ``s``,
    and ``e`` as keys, with corresponding parameter dictionaries as values.
    Each parameter dictionary includes:

    *value* float, current value of the parameter

    *step* float, step size used in the fitting procedure

    *lim* 2-list [*min*, *max*], parameter bounds

    *error* float, estimated uncertainty of the parameter

    Optionally, a parameter dictionary may also include a *tie* entry defining
    a constraint expression.

    The list :param:`exafsfit_params` includes an additional dictionary that
    defines the ``s0`` parameter (the multi-electron damping factor), along
    with optional metavariables.

    The fitting is performed using ``scipy.optimize.curve_fit()``. The
    resulting fit information returned by this function is stored in the
    :param:`exafsfit_result` attribute of the data container.
    """

    name = "EXAFS fit"
    ref = "nogui.html#exafs-fit-of-k-and-r"
    tooltip = "EXAFS fit in filtered k- and/or in r-space"
    dataAttrs = dict(x='bftk', y='bft', fit='bftfit',
                     x2='r', y2='ft', fit2='ftfit')
    plotParams = dict(fit=dict(linewidth=1.8, linestyle=':'),
                      # residue=dict(linewidth=1.5, linestyle='--')
                      )
