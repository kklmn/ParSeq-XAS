# -*- coding: utf-8 -*-
u"""
No GUI: data transformations
----------------------------

Make absorption coefficient
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In all cases, the obtained absorption coefficient is unnormalized, i.e. defined
down to an unknown multiplicative constant.

.. autoclass:: MakeTrMu
.. autoclass:: MakeFYMu
.. autoclass:: MakeTEYMu
.. autoclass:: MakeHERFD

Make EXAFS function χ(k)
~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: MakeChi

Make Fourier-transformed EXAFS function χ(r)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: MakeFT

Make back-Fourier-transformed EXAFS function χ\u0303(k)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: MakeBFT

"""

__author__ = "Konstantin Klementiev"
__date__ = "28 Nov 2023"

import sys; sys.path.append('..')  # analysis:ignore

import numpy as np
if not hasattr(np, 'trapezoid'):
    np.trapezoid = np.trapz
from numpy.polynomial import Polynomial as P

from functools import partial
from scipy.optimize import curve_fit, root
from scipy.signal import butter, sosfiltfilt
from scipy.interpolate import (CubicSpline, LSQUnivariateSpline, interp1d,
                               make_lsq_spline)
# from scipy.interpolate import BSpline
try:
    from scipy.interpolate import make_smoothing_spline
except ImportError as e:
    print('Need scipy >= 1.10.0')
    raise e

# from scipy.ndimage import uniform_filter1d
# from scipy.integrate import simps

from parseq.core import transforms as ctr
from parseq.core.logger import logger
from parseq.utils import ft as uft
from parseq.utils import math as uma
from parseq.utils.constants import eV2revA  # 2m_e(eV)/(^h(eVs)c(A/s))^2

from parseq.third_party import XAFSmass

cpus = 'half'  # can be 'all' or 'half' or a number (int)
# cpus = 4


class MakeTrMu(ctr.Transform):
    r"""
    The transmission absorption coefficient is calculated as
    :math:`µ_{tr}(E)=\log(I_0/I_{tr})`.
    """

    name = 'make tr mu'
    ref = "nogui.html#make-absorption-coefficient"
    nThreads = 1
    inArrays = ['i0', 'itr', 'eraw', 'eref']
    outArrays = ['muraw']
    defaultParams = {}

    @classmethod
    def run_main(cls, data):
        posI0 = np.where(data.i0 > 0, data.i0, np.ones_like(data.i0))
        posItr = np.where(data.itr > 0, data.itr, np.ones_like(data.itr))
        data.muraw = np.log(posI0 / posItr)
        return True


class MakeFYMu(ctr.Transform):
    name = 'make PFY mu'
    ref = "nogui.html#make-absorption-coefficient"
    nThreads = 1
    inArrays = ['i0', 'ify', 'eraw', 'eref']
    outArrays = ['muraw']
    defaultParams = {}

    @classmethod
    def run_main(cls, data):
        posI0 = np.where(data.i0 > 0, data.i0, np.ones_like(data.i0))
        posIfy = np.where(data.ify > 0, data.ify, np.ones_like(data.ify))
        data.muraw = posIfy / posI0 * posI0.max()
        return True


class MakeTEYMu(ctr.Transform):
    r"""
    These versions of :math:`µ(E)` (in fluorescence or electron yield) are
    obtained by dividing the measured signal by :math:`I_0`. In addition, the
    result is multiplied by :math:`\max(I_0)` in order to keep the physical
    meaning and units of the original signal.
    """

    name = 'make TEY mu'
    ref = "nogui.html#make-absorption-coefficient"
    nThreads = 1
    inArrays = ['i0', 'iey', 'eraw', 'eref']
    outArrays = ['muraw']
    defaultParams = {}

    @classmethod
    def run_main(cls, data):
        posI0 = np.where(data.i0 > 0, data.i0, np.ones_like(data.i0))
        posIey = np.where(data.iey > 0, data.iey, np.ones_like(data.iey))
        data.muraw = posIey / posI0 * posI0.max()
        return True


class MakeHERFD(ctr.Transform):
    r"""
    The 2D count array `xes2D` (the scan axis (DCM energy *e*) vs meridional
    detector pixel) is summed within a given band ROI horizontally (i.e. the
    columns are summed) to get a 1D count array of the length of *e*. This
    array is divided by :math:`I_0` and multiplied by :math:`\max(I_0)` to get
    the HERFD absorption coefficient in count units.
    """

    name = 'make HERFD'
    ref = "nogui.html#make-absorption-coefficient"
    nThreads = cpus
    inArrays = ['i0', 'xes2D', 'eraw', 'eref']
    outArrays = ['muraw']
    defaultParams = dict(
        cutoffNeeded=True, cutoff=20000, cutoffMaxBelow=0,
        )

    # defaultParams['roi'] = dict(
    #     kind='BandROI', name='band', use=True,
    #     begin=(300, 0), end=(300, 100), width=10)
    defaultParams['roiHERFD'] = dict(
        kind='HorizontalRangeROI', name='roi', use=True, vmin=300, vmax=400)

    @classmethod
    def run_main(cls, data):
        posI0 = np.where(data.i0 > 0, data.i0, np.ones_like(data.i0))
        dtparams = data.transformParams

        xes2Dwork = np.array(data.xes2D)
        if dtparams['cutoffNeeded']:
            cutoff = dtparams['cutoff']
            xes2Dwork[xes2Dwork > cutoff] = 0
            dtparams['cutoffMaxBelow'] = xes2Dwork.max()

        roi = dtparams['roiHERFD']
        if roi['use']:
            if roi['kind'] == 'HorizontalRangeROI':
                vmin = max(int(roi['vmin']), 0)
                vmax = int(roi['vmax'])
                posIXES = xes2Dwork[:, vmin:vmax+1].sum(axis=1)
            elif roi['kind'] == 'BandROI':
                sh = xes2Dwork.shape
                xs = np.arange(sh[1])[None, :]
                ys = data.eraw[:, None]
                m = uma.get_roi_mask(roi, xs, ys)
                masked = np.where(m, xes2Dwork, 0)
                posIXES = masked.sum(axis=1)
            else:
                raise ValueError('Unknown roi kind for {0}'.format(data.alias))
        else:
            posIXES = xes2Dwork.sum(axis=1)
        data.muraw = posIXES / posI0 * posI0.max()
        return True


class MakeChi(ctr.Transform):
    r"""
    This transformation is the largest by the number of code lines (but still
    only ~600 Python lines) and by computation load. It comprises several
    sub-steps explained below.

    Edge position :math:`E_0`
    ~~~~~~~~~~~~~~~~~~~~~~~~~

    The absorption coefficient is differentiated and optionally smoothened.
    :math:`E_0` is found within the search interval :param:`e0Where` (two
    relative fractions of the spectrum energy range) by one of the three
    methods selected by :param:`e0Method`: (0) simple derivative maximum,
    (1) the derivative maximum of interpolating cubic spline, (2, default)
    the center of mass of interpolating cubic spline.

    Energy calibration
    ~~~~~~~~~~~~~~~~~~

    The energy axis can be calibrated to match the found :math:`E_0` to a
    tabulated value. The calibration can be done by a Bragg angle offset
    (default), a lattice parameter offset or as a constant energy offset
    (:param:`eShiftKind` = 0, 1, 2). The shift can be applied by moving
    :math:`E_0` to a target value (:param:`eCalibrationMethod` = 0) or applying
    a given shift to :math:`E_0` (:param:`eCalibrationMethod` = 1). The latter
    option is useful when copying a calibration shift of :math:`E_0` found for
    a metal foil to other spectra from the same beam time measurements.

    If a calibration foil was measured *simultaneously* with the sample, two
    calibration scenarios are possible. (1) The two spectra, of the sample and
    of the foil, are loaded separately by modifying their format definitions.
    The foil spectrum is calibrated and then its :math:`E_0` shift is copied to
    the sample spectrum. (2) Only the sample spectrum is loaded and in the
    :math:`E_0` determination the derivative of the reference foil spectrum is
    used (:param:`useERefCurve` = True).

    Data rebinning
    ~~~~~~~~~~~~~~

    If the energy scan was done in a continuous way with a constant slew rate,
    the resulting spectrum is typically strongly over sampled. This means that
    several experimental points fall into one :math:`dk` interval (the EXAFS
    function will be defined on a constant :math:`dk` mesh, see below). Even
    more, as the k-space and the E-space are quadratically related, the energy
    intervals corresponding to :math:`dk` become linearly (with k) larger to
    the end of the spectrum, so more and more experimental points fall into one
    :math:`dk` interval. On converting from E-space to k-space, a kind of
    interpolation will be used, which only uses the local interpolation
    polynomial between the matching experimental points. A way of using *all*
    experimental data, called *rebinning*, consists of summing the experimental
    points belonging to one :math:`dk` interval before doing the interpolation.
    The user provides :param:`rebinRegions` dictionary that defines regions --
    pre-edge, edge, post-edge and EXAFS -- by setting their borders
    (`splitters`) and bin sizes (`deltas`). The defined bins are fed to
    ``numpy.histogram()`` to perform the actual rebinning. The number of
    original bins and the redefined bins per region are reported in
    :param:`nbinOriginal` and :param:`binDistrNew`.

    Pre-edge background
    ~~~~~~~~~~~~~~~~~~~

    The pre-edge background :math:`µ_b(E)` is constructed by polynomial
    interpolation over the region specified by :param:`preedgeWhere`. The
    polynomial law is given by :param:`preedgeExps`.

    For absorption spectra measured in transmission mode, usually a Victoreen
    polynomial :math:`aE^{-3}+bE^{-4}` or a modified Victoreen polynomial
    :math:`aE^{-3}+b` is utilized, where the coefficients are found by the
    least-squares fit from ``numpy.polynomial.Polynomial`` class.

    For absorption spectra measured in fluorescence, background subtraction is
    frequently not needed. More frequently a constant shift is sufficient
    (with only power "0"). Sometimes the spectra exhibit a net growth with
    energy, which can be approximated by a linear law (powers "0" and "1").

    Self-absorption correction
    ~~~~~~~~~~~~~~~~~~~~~~~~~~

    See the description of self-absorption correction, including its history,
    :ref:`here <sacorrection>`.

    The correction is defined by a dictionary
    :param:`selfAbsorptionCorrectionDict` that specifies the following keys:
    :param:`corrChemFormula`, :param:`corrDataTable` -- one of ("Henke",
    "BrCo", "Chantler", "Chantler total"), :param:`corrCalibEnergy` -- energy
    at which the calibration constant :param:`C` is calculated,
    :param:`corrFluoEnergy` -- the fluorescence line energy,
    :param:`corrPhiDeg, corrThetaDeg, corrTauDeg` -- the observation angles
    φ, θ and τ in degrees, :param:`corrFormula` -- 'thick' or another str,
    and :param:`corrThickness` -- the edge jump for a 'thin' case.

    If the specified material has no absorption edge within the spectrum range,
    this is reported in :param:`selfAbsorptionCorrectionDict['corrJumpStr']`.
    Otherwise, it contains a str of the tabulated edge jump.

    Post-edge background and edge normalization
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The post-edge background is needed for defining the edge normalization.

    The post-edge background is constructed by polynomial interpolation over
    the region specified by :param:`postedgeWhere`. The polynomial law is given
    by :param:`postedgeExps`. The arguments for choosing the polynomial are the
    same as for the pre-edge background.

    The obtained edge height is reported in :param:`edgeJump`.

    The post-edge background can also be used to construct a "flat" view of the
    absorption coefficient, where the post-edge part is seen horizontal. This
    can be useful for the linear combination fit and the function fit of µ(E).

    Atomic-like absorption coefficient µ₀
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The definition of EXAFS function (see below) includes µ₀ -- an artificial
    absorption coefficient that the material would have without EXAFS wiggles,
    i.e. when the central atom would be isolated without neighbors. The
    absorption coefficient of such a single atom gas is most typically not
    possible to measure, but also the electronic state of the atom in a gas
    state would differ from that in a solid or liquid, so that µ₀ has to be
    constructed artificially. µ₀ is believed to be a smooth function of energy
    and therefore is usually constructed as a spline.

    Two methods are offered for the spline creation: 'through internal k-spaced
    knots' and 'smoothing spline', :param:`mu0method` 0 and 1. The method of a
    spline through knots is generally better, as it contains only low-frequency
    oscillations, so that χ(k) preserves all the true structural oscillations.
    This is not the case for the smoothing spline, where the resulting χ(k)
    partially loses its signal and transfers it toward µ₀, but this method is
    easier to use, it always "just works".

    Prior to constructing µ₀, a helping curve is built, "µ₀ prior", that makes
    the final µ₀ resemble an absorption edge, optionally with a white line.
    µ₀ prior is first subtracted from µ before constructing a spline and then
    added again to the resulting spline.

    In the first µ₀ method ('through internal k-spaced knots'), a given number
    of knots (:param:`mu0knots`) are equidistantly placed in k-space and an LSQ
    (Least SQuared) based fitting B-spline is calculated by
    ``scipy.interpolate.LSQUnivariateSpline()``. The difference µ₀ -- µ₀ prior
    is optionally weighted with a :math:`k^w` factor (:math:`w` =
    :param:`mu0kpow`). Additionally, a given number of the first knots can be
    set variable in height to automatically minimize the low-r portion (set by
    :param:`ftMinRange`) of the FT EXAFS. The number of the varied knots
    (:param:`ftMinNKnots`) is advised to be kept small (much smaller than the
    total number of knots) to make the minimization stable.

    The second µ₀ method ('smoothing spline') depends on a smoothing parameter
    (:param:`mu0smoothingFactor`) that is set by examining the low-r FT. By
    comparing with the first µ₀ method, one can discover that the 1st FT peak
    height is always lower with the second method. This signal loss can be
    tolerated if it is smaller than the fitting error of the first shell
    coordination number.

    k-mesh and χ(k)
    ~~~~~~~~~~~~~~~

    Finally, the EXAFS function χ(k) is defined as:

    .. math::
        χ(k) = \frac{µ - µ_b - µ_0}{µ_0 - µ_b}

    where

    .. math::
        k = \sqrt{2m_e(E-E_0)} / \hbar

    The resulted χ(k) is interpolated onto the equidistant k-mesh defined by
    :param:`krange` and :param:`dk` and weighted by :math:`k^w`, where
    :math:`w` is defined by :param:`kw`.

    Denoising
    ~~~~~~~~~

    The curve can optionally be denoised by means of ``scipy.signal.butter()``.
    The two parameters, *order* and *lowpass frequency* are the first two
    parameters of ``scipy.signal.butter()``. The former one is directly it, and
    the latter is slightly modified: Wn = "lowpass frequency" × kmax. The
    calculated noise level is a normalized difference between the original and
    the denoised χ·kᵂ:

    .. math::
        N/S = \left( \sum_k(χ·k^w-χ(denoised)·k^w)^2 /
                    \sum_k(χ(denoised)·k^w)^2 \right)^{1/2}.

    """

    name = 'make chi'
    ref = "nogui.html#make-exafs-function-k"
    defaultParams = dict(
        e0Smooth=True, e0SmoothN=6, e0Where=[0.02, 0.7], e0Method=2,
        e0=None, preedgeWhere=[0.03, 0.53], preedgeExps=[-3, 0],
        postedgeWhere=[40, 400], postedgeExps=[-2, -1], edgeJump=0,
        mu0PriorIncludeWhiteLine=False, mu0PriorVScale=1., mu0PriorSmoothN=5,
        mu0method=1,  # see names in mu0methods
        mu0knots=7, mu0kpow=2,  # for mu0method=0
        ftMinimize=False, ftMinNKnots=4, ftMinRange=[0, 1],  # for mu0method=0
        mu0smoothingFactor=2e6,  # for mu0method=1
        selfAbsorptionCorrectionNeeded=False,
        selfAbsorptionCorrectionDict=dict(
            corrChemFormula='Fe2O3', corrDataTable='Chantler',
            corrCalibEnergy=7150.4, corrFluoEnergy=6404,
            corrPhiDeg=45., corrThetaDeg=45., corrTauDeg=0.,
            corrFormula='thick', corrThickness=1.),
        rebinNeeded=False,
        rebinRegions=dict(
            deltas=(1., 0.2, 0.5, 0.025),  # pre, edge, post, dk
            splitters=(-15, 15, 2.5, 'inf')),  # edge from to, k from to
        nbinOriginal=None, nbinNew=None, binDistrNew=None,
        needECalibration=False, eCalibrationMethod=1, eRef=8979.,
        eShift=0, eShiftKind=0,  # see names in eShiftKinds
        useERefCurve=False,
        kw=2, krange=[2.0, None], dk=0.025, datakmax=15.,
        denoiseNeeded=False, denoiseOrder=7, denoiseFrequency=7.0,
        noiseLevel=None,
        )
    dontSaveParamsWhenUnused = {  # paramName: paramSwitch
        'rebinRegions': 'rebinNeeded',
        'selfAbsorptionCorrectionDict': 'selfAbsorptionCorrectionNeeded'}
    nThreads = cpus  # twice as fast than nProcesses
    # nProcesses = cpus
    inArrays = ['muraw', 'eraw', 'eref']
    outArrays = ['e', 'mu', 'mu_der', 'erefrb', 'eref_der', 'e0', 'pre_edge',
                 'post_edge', 'edge_step', 'norm', 'flat', 'mu0prior', 'mu0',
                 'mu0eknots', 'k', 'chi', 'bft',
                 'mu0eknotsVaried',
                 ]
    mu0methods = ['through internal k-spaced knots', 'smoothing spline']
    eShiftKinds = ['angular shift', 'lattice shift', 'energy shift']
    firstKnot = 1

    @classmethod
    @logger(minLevel=20, attrs=[(0, 'name')])
    def get_e0(cls, data):
        dtparams = data.transformParams
        ge = np.gradient(data.e)
        good = ge > 0
        gmu = np.gradient(data.mu)
        data.mu_der = np.zeros_like(ge)
        data.mu_der[good] = gmu[good] / ge[good]
        if data.eref is not None:
            gref = np.gradient(data.erefrb)
            data.eref_der = np.zeros_like(ge)
            data.eref_der[good] = gref[good] / ge[good]
        else:
            data.eref_der = None

        if dtparams['e0Smooth']:
            ns = dtparams['e0SmoothN']
            if ns:
                data.mu_der[:] = uma.smooth_cumsum(data.mu_der, ns)
                # data.mu_der[:] = uma.smooth_savgol(data.mu_der, ns)
                data.mu_der[:ns+1] = 0.
                data.mu_der[-ns:] = 0.
                if data.eref is not None:
                    data.eref_der[:] = uma.smooth_cumsum(data.eref_der, ns)
                    # data.eref_der[:] = uma.smooth_savgol(data.eref_der, ns)
                    data.eref_der[:ns+1] = 0.
                    data.eref_der[-ns:] = 0.

        de = data.e[-1] - data.e[0]
        try:
            eminE0, emaxE0 = [data.e[0] + de*d for d in dtparams['e0Where']]
        except Exception:
            eminE0, emaxE0 = [data.e[0] + de*d
                              for d in cls.defaultParams['e0Where']]
        cond = (eminE0 <= data.e) & (data.e <= emaxE0)
        e = data.e[cond]
        # mu = data.mu[cond]
        if dtparams['useERefCurve'] and data.eref_der is not None:
            mu_der = data.eref_der[cond]
        else:
            mu_der = data.mu_der[cond]

        # e0Method=0:  simple derivative maximum
        # e0Method=1:  derivative maximum of cubic spline
        # e0Method=2:  center of mass of spline derivative
        if dtparams['e0Method'] == 0:
            e0 = e[np.argmax(mu_der)]
        elif dtparams['e0Method'] in (1, 2):
            # f = CubicSpline(e, mu)
            # der = f.derivative()
            # der2 = f.derivative(2)
            if not np.all(np.diff(e) > 0):
                raise ValueError('unsorted energy array at E={0}'.format(
                    e[1:][np.diff(e) <= 0]))
            der = CubicSpline(e, mu_der)
            der2 = der.derivative()
            xs = der2.roots()
            ys = der(xs)
            if dtparams['e0Method'] == 1:
                e0 = xs[np.argmax(ys)]
            elif dtparams['e0Method'] == 2:
                e0 = e[np.argmax(mu_der)]
                peak = np.max(ys)
                xPs = der.solve(peak*0.5)
                xleft = xPs[xPs < e0]
                xright = xPs[xPs > e0]
                if len(xleft) > 0 and len(xright) > 0:
                    x1, x2 = xleft[-1], xright[0]
                    dere = CubicSpline(e, mu_der*e)
                    e0t = dere.integrate(x1, x2) / der.integrate(x1, x2)
                    if x1 < e0t < x2:
                        e0 = e0t
        else:
            raise ValueError('unknown e0Method')
        return e0

    @classmethod
    @logger(minLevel=20, attrs=[(0, 'name')])
    def rebin(cls, data, e0):
        dtparams = data.transformParams

        rebinRegions = dtparams['rebinRegions']
        deltas = rebinRegions['deltas']
        splitters = rebinRegions['splitters']

        emin = data.e[0]  # - deltas[0]
        e1 = e0 + splitters[0]
        bins_pre = np.arange(emin, e1, deltas[0])

        e2 = e0 + splitters[1]
        bins_edge = np.arange(e1, e2, deltas[1])

        kmin = splitters[2]
        ekmin = e0 + kmin**2/eV2revA
        if e2 < ekmin:
            bins_post = np.arange(e2, ekmin, deltas[2])
        else:
            bins_post = []
            ekmin = e2
            kmin = abs((ekmin - e0)*eV2revA)**0.5

        emax = data.e[-1]
        kmaxD = abs((emax - e0)*eV2revA)**0.5
        if splitters[3] in ('inf', float('inf')):
            kmax = kmaxD
        elif isinstance(splitters[3], (float, int)):
            kmax = min(splitters[3], kmaxD)
        else:
            raise ValueError('cannot interpret kmax')
        bins_k = np.arange(kmin, kmax, deltas[3])
        bins_ke = e0 + bins_k**2/eV2revA

        bins0 = np.array([emin, e1, e2, ekmin, bins_ke[-1]])
        bins = np.array([*bins_pre, *bins_edge, *bins_post, *bins_ke])
        if not np.all(np.diff(bins) > 0):
            raise ValueError("The array of bins in not monotonic!")

        try:
            dtparams['nbinOriginal'] = np.histogram(data.e, bins0)[0]
            histNorm = np.histogram(data.e, bins)[0]
            good = histNorm > 0
            binDistrNew = []

            histNormPart = histNorm[:len(bins_pre)]
            if len(histNormPart):
                binDistrNew.append([histNormPart.min(), histNormPart.max()])
            pos = len(bins_pre)
            histNormPart = histNorm[pos:pos+len(bins_edge)]
            if len(histNormPart):
                binDistrNew.append([histNormPart.min(), histNormPart.max()])
            pos += len(bins_edge)
            histNormPart = histNorm[pos:pos+len(bins_post)]
            if len(histNormPart):
                binDistrNew.append([histNormPart.min(), histNormPart.max()])
            pos += len(bins_post)
            histNormPart = histNorm[pos:]
            if len(histNormPart):
                binDistrNew.append([histNormPart.min(), histNormPart.max()])
            dtparams['binDistrNew'] = binDistrNew
            histNormSliced = histNorm[good]
        except ValueError:
            dtparams['binDistrNew'] = None
            good = None
            histNormSliced = 1.

        # histi0 = np.histogram(data.e, bins, weights=data.i0)[0]
        # i0 = histi0[good] / histNormSliced
        # histitr = np.histogram(data.e, bins, weights=data.itr)[0]
        # itr = histitr[good] / histNormSliced
        # data.mu = np.log(i0 / itr)

        histmu = np.histogram(data.e, bins, weights=data.mu)[0]
        datamu = histmu[good] / histNormSliced

        if data.eref is not None:
            histeref = np.histogram(data.e, bins, weights=data.eref)[0]
            dataerefrb = histeref[good] / histNormSliced

        histe = np.histogram(data.e, bins, weights=data.e)[0]
        datae = histe[good] / histNormSliced

        #  change them simultaneously:
        if data.eref is not None:
            data.e, data.mu, data.erefrb = datae, datamu, dataerefrb
        else:
            data.e, data.mu = datae, datamu

        try:
            dtparams['nbinNew'] = np.histogram(data.e, bins0)[0]
        except ValueError:
            dtparams['nbinNew'] = None

    @classmethod
    @logger(minLevel=20, attrs=[(0, 'name')])
    def polyfit(cls, e, mu, exps, data):
        minPow = min(exps)
        deg = [d-minPow for d in exps]
        p = P.fit(e, mu*e**(-minPow), deg, domain=[])
        rese = p(data.e) * data.e**minPow
        rese0 = p(data.e0) * data.e0**minPow
        return rese, rese0

    @classmethod
    @logger(minLevel=20, attrs=[(0, 'name')])
    def get_pre(cls, data):
        dtparams = data.transformParams
        defpr = cls.defaultParams['preedgeWhere']
        if not dtparams['preedgeWhere']:
            dtparams['preedgeWhere'] = defpr
        emin, emax = [data.e[0] + i*(data.e0 - data.e[0]) for i in
                      dtparams['preedgeWhere']]
        if (data.e[-1] < emin) or (emax < data.e[0]):
            emin, emax = [data.e[0] + i*(data.e0 - data.e[0]) for i in defpr]
        cond = (emin <= data.e) & (data.e <= emax)
        e, mu = data.e[cond], data.mu[cond]
        return cls.polyfit(e, mu, dtparams['preedgeExps'], data)

    @classmethod
    @logger(minLevel=20, attrs=[(0, 'name')])
    def get_post(cls, data):
        dtparams = data.transformParams
        defpo = cls.defaultParams['postedgeWhere']
        if not dtparams['postedgeWhere']:
            dtparams['postedgeWhere'] = defpo
        emin, emax = [data.e0 + i for i in dtparams['postedgeWhere']]
        if (data.e[-1] < emin) or (emax < data.e[0]):
            emin, emax = [data.e0 + i for i in defpo]
        cond = (emin <= data.e) & (data.e <= emax)
        e, mu = data.e[cond], data.mu[cond]-data.pre_edge[cond]
        rese, rese0 = cls.polyfit(e, mu, dtparams['postedgeExps'], data)
        rese += data.pre_edge
        return rese, rese0

    @classmethod
    def run_main(cls, data):
        dtparams = data.transformParams

        if hasattr(data, 'eraw'):  # may be absent in data combinations
            data.e = np.array(data.eraw)
        else:
            if not hasattr(data, 'e'):  # if data combinations fails:
                return
        if hasattr(data, 'muraw'):  # may be absent in data combinations
            data.mu = np.array(data.muraw)

        if hasattr(data, 'eref'):  # may be absent in data combinations
            if data.eref is not None:
                data.erefrb = np.array(data.eref)
        else:
            data.eref = None

        # in case the analysis fails:
        data.edge_step = 1.
        data.pre_edge = np.zeros_like(data.e)

        data.e0 = cls.get_e0(data)
        if dtparams['needECalibration']:
            eRef = dtparams['eRef']
            if dtparams['eCalibrationMethod'] == 0:  # assign Eref to E0
                eShift = eRef - data.e0
                dtparams['eShift'] = eShift
            elif dtparams['eCalibrationMethod'] == 1:  # apply the shift
                eShift = dtparams['eShift']

            if dtparams['eShiftKind'] == 0:  # angular shift
                data.e = 1 / (1./data.e + 1./eRef - 1./(eRef-eShift))
            elif dtparams['eShiftKind'] == 1:  # lattice shift
                data.e *= 1 + eShift/(eRef-eShift)
            elif dtparams['eShiftKind'] == 2:  # energy shift
                data.e += eShift
            data.e0 = cls.get_e0(data)

        dtparams['nbinOriginal'] = None
        dtparams['nbinNew'] = None
        dtparams['binDistrNew'] = None
        if dtparams['rebinNeeded']:
            cls.rebin(data, data.e0)
            data.e0 = cls.get_e0(data)

        data.pre_edge, pre_e0 = cls.get_pre(data)

        res = cls.calc_self_abs(data)
        if res is not None:
            data.mu = res
            data.e0 = cls.get_e0(data)
            data.pre_edge, pre_e0 = cls.get_pre(data)

        dtparams['e0'] = data.e0

        data.post_edge, post_e0 = cls.get_post(data)

        data.edge_step = post_e0  # - pre_e0
        dtparams['edgeJump'] = data.edge_step
        data.norm = (data.mu-data.pre_edge) / data.edge_step
        data.flat = (data.mu-data.pre_edge) / (data.post_edge-data.pre_edge)

        # make mu0prior:
        data.mu0prior = np.array(data.mu)
        ie0 = np.argwhere(data.e > data.e0).flatten()[0]  # 1st point after E0
        if dtparams['mu0PriorIncludeWhiteLine']:
            ibeforeWL = np.argwhere(
                data.mu[ie0:] > data.post_edge[ie0:]).flatten()[0]
            ind = ie0 + ibeforeWL
            iafterWL = np.argwhere(
                data.mu[ind:] < data.post_edge[ind:]).flatten()[0]
            icorner = ind + iafterWL
        else:
            # build a linear rise at the edge:
            ledge = (data.e-data.e[ie0-2]) / (data.e[ie0+2]-data.e[ie0-2]) *\
                (data.mu[ie0+2]-data.mu[ie0-2]) + data.mu[ie0-2]
            try:
                icorner = np.argwhere(ledge > data.post_edge).flatten()[0]
            except IndexError:
                return
            data.mu0prior[ie0:icorner] = ledge[ie0:icorner]
        data.mu0prior[icorner:] = data.post_edge[icorner:]
        if dtparams['mu0PriorVScale'] != 1:
            data.mu0prior = \
                (data.mu0prior-data.pre_edge) * dtparams['mu0PriorVScale']\
                + data.pre_edge

        ns = dtparams['mu0PriorSmoothN']
        if ns:
            data.mu0prior[: icorner+2*ns] = \
                uma.smooth_cumsum(data.mu0prior[: icorner+2*ns], ns)
            # uma.smooth_cumsum(data.mu0prior[icorner-ns: icorner+ns], ns)
            # cornere = data.e[icorner-ns: icorner+ns+1]
            # cornermup = data.mu0prior[icorner-ns: icorner+ns+1]
            # p = P.fit(cornere, cornermup, 2, domain=[])
            # data.mu0prior[icorner-ns: icorner+ns+1] = p(cornere)

        data.mu0 = np.array(data.mu0prior)
        kmin, kmax = dtparams['krange']
        kmaxE = abs((data.e[-1] - data.e0)*eV2revA)**0.5
        kmax = min(kmax, kmaxE) if kmax else kmaxE
        dtparams['datakmax'] = kmaxE
        dk = dtparams['dk']
        data.k = np.arange(kmin, kmax + dk*0.5, dk)
        # data.k = np.arange(0, kmax + dk*0.5, dk)

        funFit = data.mu - data.mu0prior
        w = np.ones_like(data.e)  # must be positive
        wherePre = data.e < data.e0
        w[wherePre] = 1e2
        whereMax = data.e > data.e0 + kmax**2/eV2revA
        if dtparams['mu0method'] == 0:  # 'through internal k-spaced knots'
            ke = np.sign(data.e-data.e0) * (abs(data.e-data.e0)*eV2revA)**0.5
            nKnots = max(dtparams['mu0knots'], 3)
            # knots = np.linspace(kmin, kmax, nKnots)
            knots = np.linspace(0, kmax, nKnots)
            kpow = dtparams['mu0kpow']
            # above = data.e > data.e0 + kmin**2/eV2revA
            above = data.e > data.e0
            w[above] = ke[above]**kpow
            w[whereMax] = 1e-10
            try:
                spl = LSQUnivariateSpline(ke+1e-6, funFit, knots, w, ext=3)
                # spl = make_lsq_spline(ke+1e-6, funFit, knots, w=w)
            except ValueError:
                argsort = ke.argsort()
                ke = ke[argsort]
                funFit = funFit[argsort]
                spl = LSQUnivariateSpline(ke+1e-6, funFit, knots, w, ext=3)
                # spl = make_lsq_spline(ke+1e-6, funFit, knots, w=w)

            interpPrior = interp1d(ke+1e-6, data.mu0prior, assume_sorted=True)
            eknots = data.e0 + knots**2/eV2revA
            yknots = spl(knots)
            data.mu0eknots = np.array(eknots), yknots+interpPrior(knots)

            if dtparams['ftMinimize']:
                nvKnots = dtparams['ftMinNKnots']
                ftMinRange = dtparams['ftMinRange']
                r = np.fft.rfftfreq(MakeFT.nfft, dk/np.pi)
                wherer = (ftMinRange[0] <= r) & (r <= ftMinRange[1])
                # fitr = np.concatenate((r[wherer], r[wherer]))
                fitr = r[wherer]

                p0y = np.array(yknots[cls.firstKnot:nvKnots+cls.firstKnot])
                # dp = (funFit.max() - funFit.min())
                boundsy = (p0y-abs(p0y)*10, p0y+abs(p0y)*10)
                popt = curve_fit(partial(
                    cls.mu0_spline_fit, eknots=eknots, yknots=yknots,
                    e=data.e, e0=data.e0, mu0prior=data.mu0prior,
                    mu=data.mu, pre_edge=data.pre_edge, k=data.k,
                    kw=dtparams['kw'], wherer=wherer, alias=data.alias),
                    fitr, np.zeros_like(fitr), p0=p0y, bounds=boundsy)[0]
                yknots[cls.firstKnot:nvKnots+cls.firstKnot] = popt
                mu0spl = CubicSpline(eknots, yknots)
                data.mu0[above] = mu0spl(data.e[above]) + data.mu0prior[above]
                data.mu0eknotsVaried = \
                    eknots[cls.firstKnot:nvKnots+cls.firstKnot], \
                    yknots[cls.firstKnot:nvKnots+cls.firstKnot] + \
                    interpPrior(knots[cls.firstKnot:nvKnots+cls.firstKnot])
            else:
                # if True:  # both ways work equally
                data.mu0 = spl(ke) + data.mu0prior
                # else:
                #     cspl = spl.get_coeffs()
                #     # # +2 end knots + 2×order(=3) boundary knots:
                #     k = 3
                #     knotsBSpl = np.array(
                #         [ke[0]]*(k+1) + list(knots) + [ke[-1]]*(k+1))
                #     splK = BSpline(knotsBSpl, cspl, k)
                #     data.mu0 = splK(ke) + data.mu0prior
                data.mu0eknotsVaried = None

        elif dtparams['mu0method'] == 1:  # 'smoothing spline'
            data.mu0eknots = None
            data.mu0eknotsVaried = None
            s = dtparams['mu0smoothingFactor']
            if s == 0:
                s = None
            w[whereMax] = 1e-10
            try:
                spl = make_smoothing_spline(data.e, funFit, w, lam=s)
            except ValueError:
                argsort = data.e.argsort()
                data.e = data.e[argsort]
                funFit = funFit[argsort]
                spl = make_smoothing_spline(data.e, funFit, w, lam=s)
            data.mu0 = spl(data.e) + data.mu0prior
        else:
            raise ValueError(
                "unknown value mu0method={0}".format(dtparams['mu0method']))

        data.chi = cls.get_chi(
            data.e, data.e0, data.mu, data.mu0, data.pre_edge, data.k,
            dtparams['kw'], dtparams)

        # # test with ft + bft
        # # differs from VIPER by sqrt(2/pi) that is tranferred to BFT:
        # ft = np.fft.rfft(data.chi, n=MakeFT.nfft) * dk/2
        # data.bft = np.fft.irfft(ft)[0:len(data.k)] / (dk/2)

        return True

    @classmethod
    @logger(minLevel=20, attrs=[(0, 'name')])
    def calc_self_abs(cls, data):

        def fToSolve(mux):
            locExpPow = (mux+bknd)*d/corrSinPhi + sumSigma2f*d/corrSinTheta
            expFact = 1 if locExpPow.min() > 1e3 else 1 - np.exp(-locExpPow)
            return If*(mux+mubf) - corrC*mux*expFact

        dtparams = data.transformParams
        saDict = dtparams['selfAbsorptionCorrectionDict']
        if not dtparams['selfAbsorptionCorrectionNeeded']:
            saDict['corrJumpStr'] = ''
            return

        corrChemFormula = saDict['corrChemFormula']
        corrTable = saDict['corrDataTable']
        corrCalibE = saDict['corrCalibEnergy']
        corrFluoE = saDict['corrFluoEnergy']
        corrPhiDeg = saDict['corrPhiDeg']
        corrThetaDeg = saDict['corrThetaDeg']
        corrTauDeg = saDict['corrTauDeg']
        corrFormula = saDict['corrFormula']

        res = XAFSmass.parse_compound(corrChemFormula)
        if not isinstance(res, tuple):
            raise ValueError(
                "Wrong chemical formula for {0}!".format(data.alias))
        parsed, mass = res
        saDict['corrChemFormulaM'] = mass
        res = XAFSmass.calculate_element_dict(parsed, corrFluoE, corrTable)
        if isinstance(res, tuple):
            elemDictf, sumSigma2f = res[0:2]
        else:
            return
        res = XAFSmass.calculate_element_dict(parsed, corrCalibE, corrTable)
        if isinstance(res, tuple):
            elemDictc, sumSigma2c = res[0:2]
        else:
            return
        bknd = XAFSmass.calculate_absorption_background(elemDictc, data.e)

        jump, jumpElement = 0, None
        for elem, elemContr in elemDictc.items():
            if elemContr[2] > jump:
                jump = elemContr[2] * elemContr[3]  # (Δσ)*formula_coeff
                jumpElement = elem
        if jump == 0:
            saDict['corrJumpStr'] = 'no edge found'
            raise ValueError("No absorption edge for {0}!".format(data.alias))
        else:
            ss = "Δσ[{0}] (cm²/mol) = {1:.3g}".format(jumpElement, jump)
            for r in (("e-0", "e-"), ("e+0", "e+")):
                ss = ss.replace(*r)
            saDict['corrJumpStr'] = ss

        corrSinPhi = np.sin(np.radians(corrPhiDeg))
        corrSinTheta = np.sin(np.radians(corrThetaDeg))*np.cos(
            np.radians(corrTauDeg))
        if corrSinTheta < 1e-20:
            if corrSinPhi < 1e-20:
                corrG = 1
            else:
                # 'Specify non-zero escape angle'
                return
        else:
            corrG = corrSinPhi / corrSinTheta
        if (corrSinTheta < 1e-20) or (corrSinPhi < 1e-20):
            corrFormula = 'thick'
            # SendDlgItemMessage(Parent^.HWindow,id_CorrRBThick,bm_SetCheck,1,0)

        If = data.mu - data.pre_edge
        ICalib = np.interp(corrCalibE, data.e, If)
        mubf = bknd + sumSigma2f*corrG
        if corrFormula == 'thick':
            corrC = ICalib * (sumSigma2c + sumSigma2f*corrG) / jump
            muX = mubf*If / (corrC - If)
        else:
            # corrThicknessWhich = saDict['corrThicknessWhich']
            corrThickness = saDict['corrThickness']
            # if corrThicknessWhich == 0:
            #     d = corrThickness*corrSinPhi / sumSigma2c
            # elif corrThicknessWhich == 1:
            #     d = corrThickness*corrSinPhi / jump
            # else:
            #     raise ValueError("Unknown definition of sample thickeness")
            d = corrThickness*corrSinPhi / jump
            calibExpPow = sumSigma2c*d/corrSinPhi + sumSigma2f*d/corrSinTheta
            expFact = 1 - np.exp(-calibExpPow)
            corrC = ICalib * (sumSigma2c + sumSigma2f*corrG) / jump / expFact

            guess = 2 * mubf
            # calc_self_abs() has taken:
            # answ = root(fToSolve, guess, method='hybr')  # 3.178245s
            # answ = root(fToSolve, guess, method='lm')  # 9.654587s
            # answ = root(fToSolve, guess, method='broyden1')  # 0.052743s
            # answ = root(fToSolve, guess, method='broyden2')  # 0.048696s
            # answ = root(fToSolve, guess, method='anderson')  # 0.060855s
            # answ = root(fToSolve, guess, method='linearmixing')  # 0.058723s
            # answ = root(fToSolve, guess, method='diagbroyden')  # 0.058723s
            answ = root(fToSolve, guess, method='excitingmixing')  # 0.040009s
            # answ = root(fToSolve, guess, method='krylov')  # fails
            # answ = root(fToSolve, guess, method='df-sane')  # fails
            if answ.success:
                muX = answ.x
            else:
                raise ValueError(answ.message)
        return muX/jump*ICalib + data.pre_edge  # normalize it back to If

    @classmethod
    def e_to_k(cls, data, e):
        sign = np.sign(e - data.e0)
        return sign*(np.abs(e - data.e0)*eV2revA)**0.5

    @classmethod
    @logger(minLevel=20, attrs=[(0, 'name')])
    def get_chi(cls, e, e0, mu, mu0, pre_edge, k, kw, dtparams=None):
        wherek = e >= e0
        kexp = ((e[wherek] - e0)*eV2revA)**0.5
        mu_ = mu[wherek]
        mu0_ = mu0[wherek]
        pre = pre_edge[wherek]
        chie = (mu_ - mu0_) / (mu0_ - pre)
        chik = np.interp(k, kexp, chie) * k**kw
        if dtparams is not None and dtparams['denoiseNeeded']:
            sos = butter(dtparams['denoiseOrder'],
                         dtparams['denoiseFrequency']*k[-1],
                         fs=len(k), output='sos')  # Butterworth
            d = sosfiltfilt(sos, chik)
            dtparams['noiseLevel'] = (((chik-d)**2).sum() / (d**2).sum())**0.5
            chik = d
        return chik

    @classmethod
    def mu0_spline_fit(cls, r, *vals, eknots, yknots, e, e0, mu0prior, mu,
                       pre_edge, k, kw, wherer, alias):
        nvKnots = len(vals)
        yknots[cls.firstKnot:nvKnots+cls.firstKnot] = vals
        mu0spl = CubicSpline(eknots, yknots)
        mu0 = np.array(mu0prior)
        above = e > e0
        mu0[above] = mu0spl(e[above]) + mu0prior[above]
        chi = cls.get_chi(e, e0, mu, mu0, pre_edge, k, kw)  # * ftwindow
        chi -= np.trapezoid(chi, x=k) / np.trapezoid(np.ones_like(chi), x=k)
        dk = k[1] - k[0]
        ft = np.fft.rfft(chi, n=MakeFT.nfft) * dk/2
        # res = np.concatenate((ft.real[wherer]**2, ft.imag[wherer]**2))
        res = np.abs(ft[wherer])**4
        return res


class MakeFT(ctr.Transform):
    r"""
    The reason for setting χ(k) on a uniform grid is the use of the efficient
    Fast Fourier Transform (FFT) algorithm. ParSeq-XAS utilizes numpy function
    ``fft.rfft()`` that computes the one-dimensional FFT for real input.
    Because the EXAFS FT has a doubled exponential power :math:`-2ikr`, not the
    usual :math:`-ikr`, the result is multiplied by :math:`dk/2`, not the usual
    :math:`dk`; the real-space spacing :math:`dr=π/(N\cdot dk)`. The real-space
    grid is given by numpy function ``fft.rfftfreq(N, dk/π)``. The number of
    grid points :math:`N` is a class variable ``nfft`` and equals here
    2¹³ = 8192.

    Before making FT, χ(k) is multiplied by a window function that is one of
    these choices: 'none', 'box', 'linear-tapered', 'cosine-tapered',
    'Gaussian-tapered', as set by :param:`ftWindowKind`.

    Optionally, the zeroth frequency (here, distance :math:`r`) can be removed
    by nulling the first integral of χ(k); this choice is controlled by
    :param:`forceFT0`.

    The resulting FT is cut at a selected :param:`rmax` value, mainly for the
    plotting purpose.
    """

    name = 'make FT'
    ref = "nogui.html#make-fourier-transformed-exafs-function-r"
    defaultParams = dict(
        ftWindowKind='box', ftWindowProp=[1.5, 0.05],
        rmax=8.2, forceFT0=True)
    nThreads = cpus
    inArrays = ['k', 'chi']
    outArrays = ['r', 'ft', 'ftr', 'fti', 'ftwindow']

    nfft = 8192

    @classmethod
    def run_main(cls, data):
        dtparams = data.transformParams

        # # test with simps integration
        # data.r = np.arange(0.0, 10.0, 0.02)
        # kr2 = 2.0 * data.k * data.r[:, np.newaxis]
        # ft_re = simps(data.chi*np.cos(kr2), 0.5*data.k)
        # ft_im = simps(data.chi*np.sin(kr2), 0.5*data.k)
        # data.ft = (np.abs(ft_re**2 + ft_im**2))**0.5

        kind = dtparams['ftWindowKind']
        w, vmin = dtparams['ftWindowProp']
        kmin, kmax = dtparams['krange']
        # kmaxE = dtparams['datakmax']
        # kmax = min(kmax, kmaxE) if kmax else kmaxE
        data.ftwindow = uft.make_ft_window(kind, data.k, kmin, kmax, w, vmin)
        chi = np.array(data.chi) * data.ftwindow
        if dtparams['forceFT0']:
            chi -= np.trapezoid(chi, x=data.k) / np.trapezoid(
                np.ones_like(chi), x=data.k)

        dk = dtparams['dk']
        # differs from VIPER by sqrt(2/pi) that is tranferred to BFT:
        ft = np.fft.rfft(chi, n=cls.nfft) * dk/2
        r = np.fft.rfftfreq(cls.nfft, dk/np.pi)
        # wherer = slice(None) if kind == 'none' else (r <= dtparams['rmax'])
        wherer = r <= dtparams['rmax']
        data.r = r[wherer]
        data.ft = np.abs(ft)[wherer]
        data.ftr = ft.real[wherer]
        data.fti = ft.imag[wherer]
        return True


class MakeBFT(ctr.Transform):
    """
    This class applies a window function, as per :param:`bftWindowKind`,
    :param:`bftWindowRange` and :param:`bftWindowWidth`, and calculates
    Back Fourier Transform (BFT) by numpy's ``fft.irfft()``. The resulted BFT
    is cut to the k range defined in the previous steps.
    """

    name = 'make BFT'
    ref = "nogui.html#make-back-fourier-transformed-exafs-function-k"
    defaultParams = dict(
        bftWindowKind='box', bftWindowRange=[0.5, 2.5],
        bftWindowWidth=0.5)
    nThreads = cpus
    inArrays = ['k', 'r', 'ftr', 'fti', 'ftwindow']
    outArrays = ['bft', 'bftk', 'bftwindow']

    nfft = 8192

    @classmethod
    def run_main(cls, data):
        dtparams = data.transformParams
        kind = dtparams['bftWindowKind']
        if dtparams['bftWindowRange'] is not None:
            rmin, rmax = dtparams['bftWindowRange']
            w = dtparams['bftWindowWidth']
            data.bftwindow = uft.make_ft_window(kind, data.r, rmin, rmax, w)
        else:
            data.bftwindow = np.ones_like(data.r)
        dk = dtparams['dk']
        bft = np.fft.irfft((data.ftr + 1j*data.fti)*data.bftwindow, n=cls.nfft)
        if hasattr(data, 'ftwindow'):  # may not have it if loaded from file
            ftwindow = np.array(data.ftwindow)
            ftwindow[ftwindow <= 0] = 1.
        else:
            ftwindow = 1

        if not hasattr(data, 'k'):  # may not have it if loaded from file
            kmin, kmax = dtparams['krange']
            kmaxE = dtparams['datakmax']
            kmax = min(kmax, kmaxE) if kmax else kmaxE
            dk = dtparams['dk']
            data.k = np.arange(kmin, kmax + dk*0.5, dk)
        bftr = bft.real[:len(data.k)] * 2 / (dk * ftwindow)
        bftr -= np.trapezoid(bftr, x=data.k) / np.trapezoid(
            np.ones_like(bftr), x=data.k)

        kmin, kmax = dtparams['krange']
        if kmin is None:
            kmin = 0
        if kmax is None:
            kmax = 1e20
        wherek = (kmin <= data.k) & (data.k <= kmax)
        data.bft = bftr[wherek]
        data.bftk = np.array(data.k[wherek])
        return True
