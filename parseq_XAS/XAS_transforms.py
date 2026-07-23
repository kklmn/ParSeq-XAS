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
from scipy.ndimage import map_coordinates
from scipy.interpolate import (CubicSpline, LSQUnivariateSpline, interp1d,
                               make_lsq_spline)
from skimage.transform import warp
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

# cpus = 'half'  # can be 'all' or 'half' or 'quarter' a number (int)
cpus = 1


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
        healMissingFrames=True, healedMissingFrames='',
        roiHERFD=dict(kind='HorizontalRangeROI', name='roi', use=True,
                      vmin=300, vmax=400),
        dispersionCorrection=False,
        dispersionCorrectionKind=1,  # 0: edge detection, 1: elastic band
        dispersionCorrectionThreshold=0.3, dispersionCorrectionApply=False,
        roiDispersionElastic=dict(kind='BandROI', name='ela', use=True,
                                  begin=(370, 0), end=(560, 0), width=5.0))

    @classmethod
    def linear_func(cls, x, k, b):
        return k*x + b

    @classmethod
    def shear_image(cls, image, x0, k):
        def shear(xy):
            # cols are in xy[:, 0] and rows are in xy[:, 1]
            xy[:, 1] += k * (xy[:, 0] - x0)
            return xy
        return warp(image, shear, mode='edge', preserve_range=True)

    @classmethod
    def find_skew(cls, xes2D, vmin, vmax, thr, eraw):
        inBand = xes2D[:, vmin:vmax+1]
        x0 = vmin + np.argmax(inBand.sum(axis=0))
        imax = np.argmax(inBand, axis=0)
        cmax = inBand[imax, np.arange(len(imax))]
        ithr = np.argmax(inBand > thr*cmax[None, :], axis=0)
        y = ithr.astype(float)
        z = 1. / np.where(cmax > 0, cmax, 1e-20)
        z *= 0.5 / z.max()
        x = np.arange(vmin, vmax+1)
        p, _ = curve_fit(cls.linear_func, x, y, sigma=z, absolute_sigma=True)
        y0 = map_coordinates(eraw, (y,))
        skewk, skewb = p[0], p[1]
        y = map_coordinates(eraw, (skewk*x + skewb,))
        return x0, skewk, x, y0, y, z

    @classmethod
    def run_main(cls, data):
        posI0 = np.where(data.i0 > 0, data.i0, np.ones_like(data.i0))
        dtparams = data.transformParams

        data.xes2D = np.array(data.xes2Draw)
        xes2Dwork = data.xes2D
        if dtparams['cutoffNeeded']:
            cutoff = dtparams['cutoff']
            xes2Dwork[xes2Dwork > cutoff] = 0
            dtparams['cutoffMaxBelow'] = xes2Dwork.max()
        dtparams['healedMissingFrames'] = ''
        if dtparams['healMissingFrames']:
            posIXES = xes2Dwork.sum(axis=1)
            badArgs = np.argwhere(posIXES == 0).ravel()
            nbad = len(badArgs)
            if nbad > 0:
                sh = xes2Dwork.shape
                # print(data, 'detected missing frames!', badArgs)
                goodLine = np.zeros(sh[1])
                goodLines = 0
                for badArg in badArgs:
                    for goodArg in range(badArg, sh[0]):
                        if xes2Dwork[goodArg, :].sum() > 0:
                            goodLine += xes2Dwork[goodArg, :]
                            goodLines += 1
                            break
                    for goodArg in range(badArg, -1, -1):
                        if xes2Dwork[goodArg, :].sum() > 0:
                            goodLine += xes2Dwork[goodArg, :]
                            goodLines += 1
                            break
                    if goodLines > 0:
                        xes2Dwork[badArg, :] = goodLine / goodLines
                posIXES = xes2Dwork.sum(axis=1)
                # badArgs2 = np.argwhere(posIXES == 0)
                # if len(badArgs2) == 0:
                #     print(data, 'missing frames healed')
                dtparams['healedMissingFrames'] = \
                    f"healed {nbad} frame{'s' if nbad > 1 else ''}"

        roi = dtparams['roiHERFD']
        if roi['use']:
            if roi['kind'] == 'HorizontalRangeROI':
                vmin = max(int(roi['vmin']), 0) + 1
                vmax = int(roi['vmax'])
                if dtparams['dispersionCorrection']:
                    if dtparams['dispersionCorrectionKind'] == 0:  # detection
                        s = cls.find_skew(
                            xes2Dwork, vmin, vmax,
                            dtparams['dispersionCorrectionThreshold'],
                            data.eraw)
                        (x0, sk, data.skewx, data.skewy0, data.skewy,
                         data.skewz) = s
                        if dtparams['dispersionCorrectionApply']:
                            data.xes2D = cls.shear_image(xes2Dwork, x0, sk)
                            xes2Dwork = data.xes2D
                            s = cls.find_skew(
                                xes2Dwork, vmin, vmax,
                                dtparams['dispersionCorrectionThreshold'],
                                data.eraw)
                            data.skewx, data.skewy0, data.skewy, data.skewz = \
                                s[2:]
                    elif dtparams['dispersionCorrectionKind'] == 1:  # elastic
                        roiD = dtparams['roiDispersionElastic']
                        if roiD['use']:
                            x1, y1 = roiD['begin']
                            x2, y2 = roiD['end']
                            sk = (y2-y1) / (x2-x1) * len(data.eraw) /\
                                (data.eraw[-1]-data.eraw[0])
                            inBand = xes2Dwork[:, vmin:vmax+1]
                            x0 = vmin + np.argmax(inBand.sum(axis=0))
                            if dtparams['dispersionCorrectionApply']:
                                data.xes2D = cls.shear_image(xes2Dwork, x0, sk)
                                xes2Dwork = data.xes2D

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
    This transformation is the most extensive in terms of both code size
    (~600 lines of Python) and computational load. It consists of several
    sub-steps, which are described below.

    Edge position :math:`E_0`
    ~~~~~~~~~~~~~~~~~~~~~~~~~

    The absorption coefficient is differentiated and may optionally be smoothed.
    The value of :math:`E_0` is determined within the search interval defined
    by :param:`e0Where` (specified as two relative fractions of the energy
    range of the spectrum). The determination is based on one of three methods
    selected via :param:`e0Method`:

    | (0) maximum of the raw derivative,
    | (1) maximum of the derivative of an interpolating spline,
    | (2, default) center of mass of the derivative of an interpolating spline.

    Energy calibration
    ~~~~~~~~~~~~~~~~~~

    The energy axis can be calibrated so that the calculated :math:`E_0`
    matches a tabulated value. Calibration can be performed using a Bragg angle
    offset (default), a lattice parameter offset, or a constant energy offset
    (:param:`eShiftKind` = 0, 1, 2). The shift can be applied either by
    adjusting :math:`E_0` to a target value (:param:`eCalibrationMethod` = 0)
    or by applying a specified offset to :math:`E_0`
    (:param:`eCalibrationMethod` = 1). The latter approach is useful when
    transferring a calibration shift -- determined, for example, from a metal
    foil -- to other spectra measured during the same beamtime.

    If a calibration foil is measured simultaneously with the sample, two
    calibration scenarios are possible:

    1. The sample and foil spectra are loaded separately by modifying their
    format definitions. The foil spectrum is calibrated, and its resulting
    :math:`E_0` shift is then applied to the sample spectrum.

    2. Only the sample spectrum is loaded, while the derivative of the
    reference foil spectrum is used for determining :math:`E_0` (by setting
    :param:`useERefCurve` = True).

    Data rebinning
    ~~~~~~~~~~~~~~

    If the energy scan is performed in continuous scanning with a constant slew
    rate, the resulting spectrum is often strongly oversampled. In such cases,
    multiple experimental points may fall within a single :math:`dk` interval
    (the EXAFS function is defined on a uniform :math:`dk` grid, see below).
    Moreover, since k-space and energy space are quadratically related, the
    energy intervals corresponding to :math:`dk` increase linearly with k. As a
    result, progressively more experimental points are grouped within a single
    :math:`dk` interval toward the high-energy end of the spectrum.

    When converting from energy space to k-space, interpolation is applied,
    typically using only the local polynomial defined by neighboring
    experimental points. To make use of all measured data, an alternative
    approach -- *rebinning* -- can be employed. This method aggregates all
    experimental points within each :math:`dk` interval prior to interpolation.

    Rebinning is controlled by the user-defined :param:`rebinRegions`
    dictionary, which specifies regions (pre-edge, edge, post-edge, and EXAFS)
    through their boundaries (`splitters`) and bin sizes (`deltas`). The
    resulting bins are passed to ``numpy.histogram()`` to perform the rebinning.

    The number of original bins and the distribution of the new bins across
    regions are reported in :param:`nbinOriginal` and :param:`binDistrNew`,
    respectively.

    Pre-edge background
    ~~~~~~~~~~~~~~~~~~~

    The pre-edge background :math:`\mu_b(E)` is constructed using polynomial
    interpolation over the region specified by :param:`preedgeWhere`. The
    polynomial form is defined by :param:`preedgeExps`.

    For absorption spectra measured in transmission mode, a Victoreen
    polynomial :math:`aE^{-3} + bE^{-4}`, or a modified version
    :math:`aE^{-3} + b`, is typically used. The coefficients are determined via
    least-squares fitting using the ``numpy.polynomial.Polynomial`` class.

    For absorption spectra measured in fluorescence mode, explicit background
    subtraction is often unnecessary. In many cases, a constant offset (power
    "0") is sufficient. When the spectrum exhibits a gradual increase with
    energy, a linear approximation (powers "0" and "1") can be applied.

    Self-absorption correction
    ~~~~~~~~~~~~~~~~~~~~~~~~~~

    See the description of self-absorption correction, including its history,
    :ref:`here <sacorrection>`.

    The correction is defined by the dictionary
    :param:`selfAbsorptionCorrectionDict`, which includes the following keys:

    - :param:`corrChemFormula` -- chemical formula of the material,

    - :param:`corrDataTable` -- absorption data source (one of "Henke", "BrCo",
      "Chantler", or "Chantler total"),

    - :param:`corrCalibEnergy` -- energy at which the calibration constant
      :math:`C` is determined,

    - :param:`corrFluoEnergy` -- fluorescence line energy,

    - :param:`corrPhiDeg`, :param:`corrThetaDeg`, :param:`corrTauDeg` --
      observation angles φ, θ, and τ (in degrees),

    - :param:`corrFormula` -- correction model ("thick" or any other string),

    - :param:`corrThickness` -- edge jump value used in the thin sample case.

    The entry :param:`selfAbsorptionCorrectionDict['corrJumpStr']` contains a
    string representing the tabulated edge jump if it could be determined.
    If the specified material does not exhibit an absorption edge within the
    energy range of the spectrum, this is reported in the same entry.

    Post-edge background and edge normalization
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The post-edge background is required for proper edge normalization.

    It is constructed by polynomial interpolation over the region specified by
    :param:`postedgeWhere`, with the polynomial form defined by
    :param:`postedgeExps`. The same considerations for selecting the polynomial
    apply as for the pre-edge background.

    The resulting edge height is reported in :param:`edgeJump`.

    In addition, the post-edge background can be used to create a "flat"
    representation of the absorption coefficient, in which the post-edge region
    appears horizontal. This representation is particularly useful for linear
    combination fitting and function fitting of :math:`\mu(E)`.

    Atomic-like absorption coefficient µ₀
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The definition of the EXAFS function (see below) includes :math:`\mu_0` --
    an artificial absorption coefficient representing the material in the
    absence of EXAFS oscillations, i.e. as if the central atom were isolated
    without neighboring atoms. Since such a single-atom gas is typically not
    experimentally accessible -- and its electronic state would differ from
    that in condensed matter -- :math:`\mu_0` must be constructed artificially.
    It is assumed to be a smooth function of energy and is therefore usually
    approximated by a spline.

    Two methods are available for constructing the spline, controlled by
    :param:`mu0method` (0 or 1): "through internal k-spaced knots" and
    "smoothing spline". The first method is generally preferred, as it
    preserves high-frequency oscillations in :math:`\chi(k)` by restricting
    :math:`\mu_0` to low-frequency components. In contrast, the smoothing
    spline may partially absorb real EXAFS signal into :math:`\mu_0`, reducing
    the amplitude of :math:`\chi(k)`. However, it is simpler to use and
    typically more robust.

    Before constructing :math:`\mu_0`, an auxiliary curve "µ₀ prior" is
    generated to shape the result to resemble an absorption edge, optionally
    including a white line. This prior is subtracted from :math:`\mu` before
    spline fitting and added back afterward.

    In the first method (spline through k-spaced knots), a specified number of
    knots (:param:`mu0knots`) are placed uniformly in k-space, and a
    least-squares B-spline is computed using
    ``scipy.interpolate.LSQUnivariateSpline()``. The difference
    :math:`\mu_0 - \mu_0^{\mathrm{prior}}` can optionally be weighted by
    :math:`k^w` (with :math:`w` = :param:`mu0kpow`). Additionally, a selected
    number of initial knots can be allowed to vary in height to minimize the
    low-R region (defined by :param:`ftMinRange`) of the Fourier-transformed
    EXAFS signal. The number of such variable knots (:param:`ftMinNKnots`)
    should be kept small relative to the total number of knots to ensure stable
    minimization.

    The second method (smoothing spline) depends on a smoothing parameter
    (:param:`mu0smoothingFactor`), which is typically adjusted by examining the
    low-R region of the Fourier transform. Compared to the first method, the
    resulting first Fourier peak is consistently lower due to partial signal
    loss. This attenuation is acceptable if it remains smaller than the
    uncertainty in determining the first-shell coordination number.

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

    The EXAFS curve can optionally be denoised using ``scipy.signal.butter()``.
    The two main parameters -- order and low-pass frequency -- correspond to
    the first two arguments of ``scipy.signal.butter()``. The former is used
    directly, while the latter is scaled as
    :math:`W_n = (\text{low-pass frequency}) \times k_{\max}`. The noise level
    is estimated as the normalized difference between the original and the
    denoised :math:`\chi \cdot k^w` curve:

    .. math::
        N/S = \left( \sum_k(χ·k^w-χ_{\rm denoised}·k^w)^2 /
                    \sum_k(χ_{\rm denoised}·k^w)^2 \right)^{1/2}.

    """

    name = 'make chi'
    ref = "nogui.html#make-exafs-function-k"
    defaultParams = dict(
        e0Smooth=True, e0SmoothN=6, e0Where=[0.02, 0.7], e0Method=2, e0=None,
        preedgeWhere=[0.03, 0.53], preedgeExps=[-3, 0], preedgePinPoint=0,
        postedgeWhere=[40, 400], postedgeExps=[-2, -1], postedgePinPoint=0,
        edgeJump=0,
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
        pinholeCorrectionNeeded=False,
        pinholeCorrectionDict=dict(fraction=0.5),
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
        'selfAbsorptionCorrectionDict': 'selfAbsorptionCorrectionNeeded',
        'pinholeCorrectionDict': 'pinholeCorrectionNeeded'}
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
        gmu = np.gradient(data.mu)
        data.mu_der = np.zeros_like(ge)
        good = ge > 0
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
        if len(mu_der) == 0:
            return data.e[0]

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
    def polyfit(cls, e, mu, exps, pinpoint, data):
        minPow = min(exps)
        deg = [d-minPow for d in exps]
        if isinstance(pinpoint, (list, tuple)) and (len(pinpoint) == 2):
            e = np.append(e, pinpoint[0])
            mu = np.append(mu, pinpoint[1])
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
        pre_emin, pre_emax = [data.e[0] + i*(data.e0 - data.e[0])
                              for i in dtparams['preedgeWhere']]
        if (data.e[-1] < pre_emin) or (pre_emax < data.e[0]):
            pre_emin, pre_emax = [data.e[0] + i*(data.e0 - data.e[0])
                                  for i in defpr]
        cond = (pre_emin <= data.e) & (data.e <= pre_emax)
        e, mu = data.e[cond], data.mu[cond]
        if len(e) == 1:
            return mu[0]*np.ones_like(data.mu), data.e0
        else:
            return cls.polyfit(e, mu, dtparams['preedgeExps'],
                               dtparams['preedgePinPoint'], data)

    @classmethod
    @logger(minLevel=20, attrs=[(0, 'name')])
    def get_post(cls, data):
        dtparams = data.transformParams
        defpo = cls.defaultParams['postedgeWhere']
        if not dtparams['postedgeWhere']:
            dtparams['postedgeWhere'] = defpo
        post_emin, post_emax = [data.e0 + i for i in dtparams['postedgeWhere']]
        if (data.e[-1] < post_emin) or (post_emax < data.e[0]):
            post_emin, post_emax = [data.e0 + i for i in defpo]
            dtparams['postedgeWhere'] = defpo
        if (data.e[-1] < post_emin) or (post_emax < data.e[0]):
            post_emin = data.e[-5]
            dtparams['postedgeWhere'][0] = post_emin - data.e0
        cond = (post_emin <= data.e) & (data.e <= post_emax)
        e, mu = data.e[cond], data.mu[cond]-data.pre_edge[cond]
        rese, rese0 = cls.polyfit(e, mu, dtparams['postedgeExps'],
                                  dtparams['postedgePinPoint'], data)
        rese += data.pre_edge
        return rese, rese0

    @classmethod
    @logger(minLevel=20, attrs=[(0, 'name')])
    def get_mu0prior(cls, data):
        dtparams = data.transformParams
        data.mu0prior = np.array(data.mu)
        ie0s = np.argwhere(data.e > data.e0).flatten()
        if len(ie0s) == 0:
            return
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

    @classmethod
    @logger(minLevel=20, attrs=[(0, 'name')])
    def get_mu0(cls, data):
        dtparams = data.transformParams
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
                r = np.fft.rfftfreq(uft.nfft, dk/np.pi)
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

        res = cls.calc_pinhole(data)
        if res is not None:
            data.mu = res
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

        cls.get_mu0prior(data)
        cls.get_mu0(data)

        data.chi = cls.get_chi(
            data.e, data.e0, data.mu, data.mu0, data.pre_edge, data.k,
            dtparams['kw'], dtparams)

        # # test with ft + bft
        # # differs from VIPER by sqrt(2/pi) that is tranferred to BFT:
        # ft = np.fft.rfft(data.chi, n=uft.nfft) * dk/2
        # data.bft = np.fft.irfft(ft)[0:len(data.k)] / (dk/2)

        return True

    @classmethod
    @logger(minLevel=20, attrs=[(0, 'name')])
    def calc_pinhole(cls, data):
        dtparams = data.transformParams
        if not dtparams['pinholeCorrectionNeeded']:
            return
        phDict = dtparams['pinholeCorrectionDict']
        x = phDict['fraction']
        muTrue = np.log(np.abs((1-x) / (np.exp(-data.mu)-x)))
        return muTrue

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
                break
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
        ft = np.fft.rfft(chi, n=uft.nfft) * dk/2
        # res = np.concatenate((ft.real[wherer]**2, ft.imag[wherer]**2))
        res = np.abs(ft[wherer])**4
        return res


class MakeFT(ctr.Transform):
    r"""
    The use of a uniform :math:`\chi(k)` grid enables efficient computation via
    the Fast Fourier Transform (FFT). ParSeq-XAS uses the NumPy function
    ``fft.rfft()`` to compute the one-dimensional FFT for real-valued input.

    Because the EXAFS Fourier transform uses the kernel :math:`\exp(-2ikr)`
    rather than the standard :math:`\exp(-ikr)`, the result is scaled by
    :math:`dk/2` instead of the usual :math:`dk`. The real-space sampling
    interval is given by :math:`dr = \pi / (N \cdot dk)`. The corresponding
    real-space grid is obtained using ``fft.rfftfreq(N, dk/π)``, where the
    number of grid points :math:`N` is defined by the class variable ``nfft``
    (set to 4096).

    Before performing the Fourier transform, :math:`\chi(k)` is multiplied by a
    window function selected via :param:`ftWindowKind`. Available options
    include 'none', 'box', 'linear-tapered', 'cosine-tapered', and
    'Gaussian-tapered'.

    Optionally, the zero-frequency component (corresponding here
    :math:`r = 0`) can be suppressed by enforcing a zero integral of
    :math:`\chi(k)k^w`. This behavior is controlled by :param:`forceFT0`.

    For visualization purposes, the resulting Fourier transform is truncated at
    a user-defined maximum distance :param:`rmax`.
    """

    name = 'make FT'
    ref = "nogui.html#make-fourier-transformed-exafs-function-r"
    defaultParams = dict(
        ftWindowKind='box', ftWindowProp=[1.5, 0.05],
        rmax=8.2, forceFT0=True)
    nThreads = cpus
    inArrays = ['k', 'chi']
    outArrays = ['r', 'ft', 'ftr', 'fti', 'ftwindow']

    nfft = uft.nfft

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
        if kmax is None:
            kmax = dtparams['datakmax']
        data.ftwindow = uft.make_ft_window(kind, data.k, kmin, kmax, w, vmin)
        chi = np.array(data.chi) * data.ftwindow
        if dtparams['forceFT0']:
            norm = np.trapezoid(np.ones_like(chi), x=data.k)
            if norm == 0:
                norm = 1
            chi -= np.trapezoid(chi, x=data.k) / norm

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
    This class applies a window function, as defined by :param:`bftWindowKind`,
    :param:`bftWindowRange`, and :param:`bftWindowWidth`, and computes the Back
    Fourier Transform (BFT) using NumPy's ``fft.irfft()``.

    The resulting BFT is then restricted to the k-range defined in the
    preceding steps.
        """

    name = 'make BFT'
    ref = "nogui.html#make-back-fourier-transformed-exafs-function-k"
    defaultParams = dict(
        bftWindowKind='box', bftWindowRange=[0.5, 2.5],
        bftWindowWidth=0.5)
    nThreads = cpus
    inArrays = ['k', 'r', 'ftr', 'fti', 'ftwindow']
    outArrays = ['bft', 'bftk', 'bftwindow']

    nfft = uft.nfft

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
        norm = np.trapezoid(np.ones_like(bftr), x=data.k)
        if norm == 0:
            norm = 1
        bftr -= np.trapezoid(bftr, x=data.k) / norm

        kmin, kmax = dtparams['krange']
        if kmin is None:
            kmin = 0
        if kmax is None:
            kmax = 1e20
        wherek = (kmin <= data.k) & (data.k <= kmax)
        data.bft = bftr[wherek]
        data.bftk = np.array(data.k[wherek])
        return True
