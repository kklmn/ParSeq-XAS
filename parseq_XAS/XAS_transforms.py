# -*- coding: utf-8 -*-
__author__ = "Konstantin Klementiev"
__date__ = "28 Nov 2023"

import sys; sys.path.append('..')  # analysis:ignore

import numpy as np
if not hasattr(np, 'trapezoid'):
    np.trapezoid = np.trapz
from numpy.polynomial import Polynomial as P

# from functools import partial
from scipy import optimize
from scipy.interpolate import CubicSpline, LSQUnivariateSpline, interp1d
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
from parseq.utils.constants import eV2revA

from parseq.third_party import XAFSmass

cpus = 'half'  # can be 'all' or 'half' or a number (int)
# cpus = 4


class MakeTrMu(ctr.Transform):
    name = 'make tr mu'
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
    name = 'make TEY mu'
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
    name = 'make HERFD'
    nThreads = cpus
    inArrays = ['i0', 'xes2D', 'eraw', 'eref']
    outArrays = ['muraw']
    defaultParams = dict(
        cutoffNeeded=True, cutoff=2000, cutoffMaxBelow=0,
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
                xs = data.eraw
                ys = np.arange(sh[0])[:, None]
                m = uma.get_roi_mask(roi, xs, ys)
                _, indx = np.nonzero(m)
                posIXES = xes2Dwork[m].sum(axis=1) if len(indx) > 0 else \
                    xes2Dwork.sum(axis=1)
            else:
                raise ValueError('Unknown roi kind for {0}'.format(data.alias))
        else:
            posIXES = xes2Dwork.sum(axis=1)
        data.muraw = posIXES / posI0 * posI0.max()
        return True


class MakeChi(ctr.Transform):
    name = 'make chi'
    defaultParams = dict(
        e0Smooth=True, e0SmoothN=6, e0Where=[0.02, 0.7], e0Method=2,
        e0=None, preedgeWhere=[0.03, 0.33], preedgeExps=[-3, 0],
        postedgeWhere=[40, 400], postedgeExps=[-2, -1], edgeJump=0,
        mu0PriorIncludeWhiteLine=False, mu0PriorVScale=1., mu0PriorSmoothN=5,
        mu0method=1,  # see names in mu0methods
        mu0knots=7, mu0kpow=2,  # for mu0method=0
        # ftMinimize=False, ftMinNKnots=4, ftMinRange=[0, 1],  # for method=0
        mu0smoothingFactor=2e6,  # for mu0method=1
        selfAbsorptionCorrectionNeeded=False,
        selfAbsorptionCorrectionDict=dict(
            corrChemFormula='Fe2O3', corrDataTable='Chantler',
            corrCalibEnergy=7150.4, corrFluoEnergy=6404,
            corrPhiDeg=45., corrThetaDeg=45., corrTauDeg=0.,
            corrFormula='thick', corrThickness=1.),
        rebinNeeded=False,
        rebinRegions=dict(
            deltas=(1., 0.2, 0.5, 0.025), splitters=(-15, 15, 2.5, 'inf')),
        nbinOriginal=None, nbinNew=None, binDistrNew=None,
        needECalibration=False, eCalibrationMethod=1, eRef=8979.,
        eShift=0, eShiftKind=0,  # see names in eShiftKinds
        useERefCurve=False,
        kw=2, krange=[2.0, None], dk=0.025, datakmax=15.
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
                 # 'mu0eknotsVaried',
                 ]
    mu0methods = ['through internal k-spaced knots', 'smoothing spline']
    eShiftKinds = ['angular shift', 'lattice shift', 'energy shift']

    RED_TXT = '<span style="font-weight:600;color:#aa0000;">{0}</span>'
    GREEN_TXT = '<span style="color:#00aa00;">{0}</span>'

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

        # e0Method = 0:  simple derivative maximum
        # e0Method = 1:  derivative maximum of cubic spline
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

        data.e = np.array(data.eraw)
        data.mu = np.array(data.muraw)
        if data.eref is not None:
            data.erefrb = np.array(data.eref)
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
        # icorner = np.argwhere(data.mu > data.post_edge).flatten()[0]
        ie0 = np.argwhere(data.e > data.e0).flatten()[0]
        # build a linear rise at the edge:
        ledge = (data.e - data.e[ie0-2]) / (data.e[ie0+2] - data.e[ie0-2]) *\
            (data.mu[ie0+2] - data.mu[ie0-2]) + data.mu[ie0-2]
        icorner = np.argwhere(ledge > data.post_edge).flatten()[0]
        if dtparams['mu0PriorIncludeWhiteLine']:
            icornerW = np.argwhere(
                data.mu[icorner:] < data.post_edge[icorner:]).flatten()[0]
            icorner += icornerW
        else:
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
            # whereMin = data.e > data.e0 + kmin**2/eV2revA
            whereMin = data.e > data.e0
            w[whereMin] = ke[whereMin]**kpow
            w[whereMax] = 1e-10
            try:
                spl = LSQUnivariateSpline(ke+1e-6, funFit, knots, w, ext=3)
            except ValueError:
                argsort = ke.argsort()
                ke = ke[argsort]
                funFit = funFit[argsort]
                spl = LSQUnivariateSpline(ke+1e-6, funFit, knots, w, ext=3)

            # cspl = spl.get_coeffs()
            interpPrior = interp1d(ke+1e-6, data.mu0prior, assume_sorted=True)
            eknots = data.e0 + knots**2/eV2revA
            data.mu0eknots = eknots, spl(knots) + interpPrior(knots)

            # !! curve_fit to minimize low-r FT is unstable !! removed
            # # +2 end knots + 2×order(=3) boundary knots:
            # knotsBSpl = np.array(4*[ke[0]] + list(knots) + 4*[ke[-1]])
            # if dtparams['ftMinimize']:
            #     ftMinNKnots = dtparams['ftMinNKnots']
            #     # 2 first intervals are fixed:
            #     variedCoeffs = cspl[2:ftMinNKnots+2]
            #     ftMinRange = dtparams['ftMinRange']
            #     r = np.fft.rfftfreq(MakeFT.nfft, dk/np.pi)
            #     wherer = (ftMinRange[0] <= r) & (r <= ftMinRange[1])
            #     # fitr = np.concatenate((r[wherer], r[wherer]))
            #     fitr = r[wherer]
            #     popt = curve_fit(partial(
            #         cls.mu0_BSpline, knots=knotsBSpl, allCoeffs=cspl,
            #         ke=ke, mu0prior=data.mu0prior, e=data.e, e0=data.e0,
            #         mu=data.mu, pre_edge=data.pre_edge, k=data.k,
            #         kw=dtparams['kw'], wherer=wherer),
            #         fitr, np.zeros_like(fitr), p0=variedCoeffs)[0]
            #     cspl[2:2+ftMinNKnots] = popt
            #     splK = BSpline(knotsBSpl, cspl, 3)
            #     data.mu0 = splK(ke) + data.mu0prior
            #     knotsK = knots[1:-1][:ftMinNKnots]
            #     knotsV = splK(knotsK) + interpPrior(knotsK)
            #     knotsE = eknots[1:-1][:ftMinNKnots]
            #     data.mu0eknotsVaried = (knotsE, knotsV)
            # else:
            #     if True:  # both ways work equally
            #         data.mu0 = spl(ke) + data.mu0prior
            #     else:
            #         splK = BSpline(knotsBSpl, cspl, 3)
            #         data.mu0 = splK(ke) + data.mu0prior
            #     data.mu0eknotsVaried = None
            data.mu0 = spl(ke) + data.mu0prior
        elif dtparams['mu0method'] == 1:  # 'smoothing spline'
            data.mu0eknots = None
            # data.mu0eknotsVaried = None
            s = dtparams['mu0smoothingFactor']
            if s == 0:
                s = None
            whereMax = data.e > data.e0 + kmax**2/eV2revA
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

        data.chi = cls.get_chi(data.e, data.e0, data.mu, data.mu0,
                               data.pre_edge, data.k, dtparams['kw'])

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
        if not dtparams['selfAbsorptionCorrectionNeeded']:
            return

        saDict = dtparams['selfAbsorptionCorrectionDict']
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
            saDict['corrJumpStr'] = cls.RED_TXT.format('no edge found')
            raise ValueError("No absorption edge for {0}!".format(data.alias))
        else:
            ss = "Δσ[{0}] (cm²/mol) = {1:.3g}".format(jumpElement, jump)
            for r in (("e-0", "e-"), ("e+0", "e+")):
                ss = ss.replace(*r)
            saDict['corrJumpStr'] = cls.GREEN_TXT.format(ss)

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
            root = optimize.root(fToSolve, guess, method='krylov')
            if root.success:
                muX = root.x
            else:
                raise ValueError(root.message)
        return muX/jump*ICalib + data.pre_edge  # normalize it back to If

    @classmethod
    @logger(minLevel=20, attrs=[(0, 'name')])
    def get_chi(cls, e, e0, mu, mu0, pre_edge, k, kw):
        wherek = e >= e0
        kexp = ((e[wherek] - e0)*eV2revA)**0.5
        mu_ = mu[wherek]
        mu0_ = mu0[wherek]
        pre = pre_edge[wherek]
        chie = (mu_ - mu0_) / (mu0_ - pre)
        return np.interp(k, kexp, chie) * k**kw

    # @classmethod
    # def mu0_BSpline(cls, r, *vals, knots, allCoeffs, ke, mu0prior, e, e0, mu,
    #                 pre_edge, k, kw, wherer):
    #     c = list(allCoeffs)
    #     c[2:2+len(vals)] = vals
    #     splK = BSpline(knots, c, 3)
    #     mu0 = splK(ke) + mu0prior
    #     chi = cls.get_chi(e, e0, mu, mu0, pre_edge, k, kw) # * ftwindow
    #     chi -= np.trapezoid(chi, x=k) / np.trapezoid(np.ones_like(chi), x=k)
    #     dk = k[1] - k[0]
    #     ft = np.fft.rfft(chi, n=MakeFT.nfft) * dk/2
    #     # return np.concatenate((ft.real[wherer], ft.imag[wherer]))
    #     return np.abs(ft[wherer])


class MakeFT(ctr.Transform):
    name = 'make FT'
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
    name = 'make BFT'
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
        ftwindow = np.array(data.ftwindow)
        ftwindow[ftwindow <= 0] = 1.
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
