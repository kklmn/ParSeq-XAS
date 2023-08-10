# -*- coding: utf-8 -*-
__author__ = "Konstantin Klementiev"
__date__ = "6 Mar 2023"

import numpy as np

from silx.gui import qt, icons
from scipy.interpolate import interp1d
from functools import partial

import sys; sys.path.append('..')  # analysis:ignore
from parseq.core import singletons as csi
from parseq.gui.dataRebin import DataRebinWidget
from parseq.gui.roi import RangeWidget, RangeWidgetSplit, RoiWidget
import parseq.gui.gcommons as gco

# from parseq.core import commons as cco
from parseq.gui.propWidget import PropWidget
from parseq.utils import ft as uft
from parseq.third_party import XAFSmass

from . import XAS_transforms as xtr


class RangeWidgetE0(RangeWidget):
    'fractionMin, fractionMax of range [eMin, eMax]'

    def convert(self, direction, vmin, vmax):
        if len(csi.selectedItems) == 0:
            return
        data = csi.selectedItems[0]
        if not hasattr(data, 'e'):
            return
        dE = data.e[-1] - data.e[0]
        if direction == 1:  # from vidget values to roi limits
            amin, amax = [data.e[0] + dE*v for v in (vmin, vmax)]
        elif direction == -1:  # from roi limits to vidget
            amin, amax = [(v - data.e[0]) / dE for v in (vmin, vmax)]
            if not (0 <= amin <= 1):
                amin = 0
            if not (0 <= amax <= 1):
                amax = 1
            if amin > amax:
                amin, amax = amax, amin
        return amin, amax


class RangeWidgetPre(RangeWidget):
    'fractionMin, fractionMax of range [eMin, E0]'

    def convert(self, direction, vmin, vmax):
        if len(csi.selectedItems) == 0:
            return
        data = csi.selectedItems[0]
        if not (hasattr(data, 'e') and hasattr(data, 'e0')):
            return
        dE0 = data.e0 - data.e[0]
        if direction == 1:  # from vidget values to roi limits
            amin, amax = [data.e[0] + dE0*v for v in (vmin, vmax)]
        elif direction == -1:  # from roi limits to vidget
            amin, amax = [(v - data.e[0]) / dE0 for v in (vmin, vmax)]
            if not (0 <= amin <= 1):
                amin = 0
            if not (0 <= amax <= 1):
                amax = 1
            if amin > amax:
                amin, amax = amax, amin
        return amin, amax


class RangeWidgetPost(RangeWidget):
    'min+E0, max+E0 (eV)', 'post-edge'

    def convert(self, direction, vmin, vmax):
        if len(csi.selectedItems) == 0:
            return
        data = csi.selectedItems[0]
        if not (hasattr(data, 'e') and hasattr(data, 'e0')):
            return
        if direction == 1:  # from vidget values to roi limits
            amin, amax = [data.e0 + v for v in (vmin, vmax)]
        elif direction == -1:  # from roi limits to vidget
            amin, amax = [v - data.e0 for v in (vmin, vmax)]
            if amin < 0:
                amin = 0
            if amin > amax:
                amin, amax = amax, amin
        return amin, amax


class RangeWidgetFTWidthAndMin(RangeWidgetSplit):
    plotFactor = 0.98

    def convert(self, direction, vmin, vmax):
        if len(csi.selectedItems) == 0:
            return
        data = csi.selectedItems[0]
        dtparams = data.transformParams
        kmin, kmax = dtparams['krange']
        if not kmax:
            kmax = dtparams['datakmax']
        ymax = self.plot.getYAxis().getLimits()[1] * self.plotFactor
        if direction == 1:  # from spinbox values to roi limits
            if vmin > (kmax - kmin)*0.5:
                vmin = (kmax - kmin)*0.5
                self.minSpinBox.setValue(vmin)
            amin = vmin + kmin if vmin else kmin
            amax = vmax*ymax if vmax else ymax
        elif direction == -1:  # from roi limits to spinbox values
            needMove = False

            amin = vmin - kmin if vmin else 0
            if amin < self.minSpinBox.minimum():
                amin = self.minSpinBox.minimum()
                needMove = True
            if amin > (kmax - kmin)*0.5:
                amin = (kmax - kmin)*0.5
                needMove = True
            if amin > self.minSpinBox.maximum():
                amin = self.minSpinBox.maximum()
                needMove = True

            amax = vmax/ymax if vmax else 1/ymax
            if amax < self.maxSpinBox.minimum():
                amax = self.maxSpinBox.minimum()
                needMove = True
            if amax > self.maxSpinBox.maximum():
                amax = self.maxSpinBox.maximum()
                needMove = True

            if needMove:
                rmin, rmax = self.convert(1, amin, amax)
                self.roi.setPosition((rmin, rmax))

        return amin, amax


class HERFDWidget(PropWidget):
    u"""
    Help page under construction

    .. image:: _images/mickey-rtfm.gif
        :width: 309

    test link: `MAX IV Laboratory <https://www.maxiv.lu.se/>`_

    """

    name = 'extract HERFD'

    def __init__(self, parent=None, node=None):
        super().__init__(parent, node)
        plot = self.node.widget.plot

        layout = qt.QVBoxLayout()

        roiPanel = qt.QGroupBox(self)
        roiPanel.setFlat(False)
        roiPanel.setTitle(u'HERFD ROI')
        layoutR = qt.QVBoxLayout()
        layoutR.setContentsMargins(0, 2, 2, 2)
        self.roiWidget = RoiWidget(
            self, plot, ['HorizontalRangeROI', 'BandROI'],
            fmt=['range: {0[0]:.0f}, {0[1]:.0f}', 'auto'])
        self.roiWidget.acceptButton.clicked.connect(self.acceptBand)
        self.registerPropWidget(
            [self.roiWidget.table, self.roiWidget.acceptButton],
            'HERFD roi', 'roiHERFD', transformNames='make HERFD')
        layoutR.addWidget(self.roiWidget)
        roiPanel.setLayout(layoutR)
        roiPanel.setSizePolicy(qt.QSizePolicy.Minimum, qt.QSizePolicy.Fixed)
        layout.addWidget(roiPanel)

        layout.addStretch()
        self.setLayout(layout)

    def acceptBand(self):
        self.roiWidget.syncRoi()
        curRoi = self.roiWidget.getCurrentRoi()
        self.updateProp('roiHERFD', curRoi)
        for data in csi.selectedItems:
            data.transformParams['roiHERFD']['use'] = True
        # nextWidget = csi.nodes[u'µd'].widget.transformWidget
        # nextWidget.setUIFromData()

    def extraSetUIFromData(self):
        if len(csi.selectedItems) == 0:
            return
        data = csi.selectedItems[0]
        if hasattr(data, 'xes2D'):
            self.roiWidget.dataToCount = data.xes2D  # to display roi counts
            dtparams = data.transformParams
            self.roiWidget.setRois(dict(dtparams['roiHERFD']))


class MuWidget(PropWidget):
    u"""
    Help page under construction

    .. image:: _images/mickey-rtfm.gif
        :width: 309

    test link: `MAX IV Laboratory <https://www.maxiv.lu.se/>`_

    """

    name = u'subtract bknd and make µ\u2080'

    properties = {
        'subtract_preedge': True, 'normalize': True, 'flatten': False,
        'show_derivative': False, 'show_E0': True,
        'show_mu0': True, 'show_mu0_prior': False, 'show_post': False}

    extraLines = 'pre_edge', 'post_edge', 'mu0', 'mu0prior', 'e0', 'mu_der'
    preEdgePlotParams = {'linewidth': 1.0, 'linestyle': '-.'}
    postEdgePlotParams = {'linewidth': 1.0, 'linestyle': '--'}
    mu0PlotParams = {'linewidth': 0.75, 'linestyle': '--'}
    mu0priorPlotParams = {'linewidth': 0.75, 'linestyle': ':'}
    derPlotParams = {'linewidth': 1.0, 'linestyle': '-', 'yaxis': 'right'}
    e0PlotParams = {'linewidth': 0.5, 'linestyle': '-'}

    cursorLabels = ['E', 'k', u'µd']

    def __init__(self, parent=None, node=None):
        super().__init__(parent, node)
        plot = self.node.widget.plot

        layout = qt.QVBoxLayout()

        layoutDM = qt.QHBoxLayout()
        self.checkBoxDerivative = qt.QCheckBox("plot µd'")
        self.checkBoxDerivative.setChecked(self.properties['show_derivative'])
        self.checkBoxDerivative.toggled.connect(
            partial(self.showSlot, 'show_derivative'))
        layoutDM.addWidget(self.checkBoxDerivative)
        self.checkBoxDerivativeRef = qt.QCheckBox("use (Eref)' instead")
        self.registerPropWidget(self.checkBoxDerivativeRef, 'use Eref',
                                'useERefCurve')
        layoutDM.addWidget(self.checkBoxDerivativeRef)
        layoutDM.addStretch()
        layout.addLayout(layoutDM)

        e0Panel = qt.QGroupBox(self)
        e0Panel.setFlat(False)
        e0Panel.setStyleSheet('QGroupBox[flat="false"] {font-weight: bold;}')
        e0Panel.setTitle('E\u2080')
        layoutE = qt.QVBoxLayout()
        layoutE.setContentsMargins(10, 2, 2, 2)

        smoothPanel = qt.QGroupBox(self)
        smoothPanel.setFlat(True)
        smoothPanel.setTitle('n-point smoothing of derivative')
        smoothPanel.setCheckable(True)
        self.registerPropWidget(smoothPanel, smoothPanel.title(), 'e0Smooth')
        layoutSm = qt.QHBoxLayout()
        layoutSm.setContentsMargins(10, 2,  2, 2)
        labeNSm = qt.QLabel('n =')
        layoutSm.addWidget(labeNSm)
        self.e0SmoothN = qt.QSpinBox()
        self.e0SmoothN.setToolTip(u'n ≥ 2')
        self.e0SmoothN.setMinimum(1)
        self.e0SmoothN.setFixedWidth(44)
        layoutSm.addWidget(self.e0SmoothN)
        layoutSm.addStretch()
        smoothPanel.setLayout(layoutSm)
        layoutE.addWidget(smoothPanel)
        self.registerPropWidget((labeNSm, self.e0SmoothN), 'n', 'e0SmoothN')

        e0WherePanel = qt.QGroupBox(self)
        e0WherePanel.setFlat(False)
        e0WherePanel.setTitle('position of E\u2080')
        layoutW = gco.QVBoxLayoutAbove()
        layoutW.setContentsMargins(10, 2, 2, 2)
        layoutW.setSpacing(2)

        checkBoxShowE0 = qt.QCheckBox('plot E\u2080', e0WherePanel)
        checkBoxShowE0.setChecked(self.properties['show_E0'])
        checkBoxShowE0.toggled.connect(partial(self.showSlot, 'show_E0'))
        layoutW.addExtraWidget(checkBoxShowE0)

        self.e0RangeWidget = RangeWidgetE0(
            self, plot, 'search range',
            'fractionMin, fractionMax of range [eMin, eMax]', 'E0-range',
            "#da70d6", "{0[0]:.3f}, {0[1]:.3f}",
            xtr.MakeChi.defaultParams['e0Where'])
        self.registerPropWidget(self.e0RangeWidget, 'E0 search range',
                                'e0Where')
        layoutW.addWidget(self.e0RangeWidget)

        e0Where0 = qt.QRadioButton(
            'simple maximum of derivative', e0WherePanel)
        layoutW.addWidget(e0Where0)
        e0Where1 = qt.QRadioButton(
            'maximum of cubic spline derivative', e0WherePanel)
        layoutW.addWidget(e0Where1)
        e0Where2 = qt.QRadioButton(
            'center of mass of spline derivative', e0WherePanel)
        layoutW.addWidget(e0Where2)
        e0WherePanel.setLayout(layoutW)
        layoutE.addWidget(e0WherePanel)
        self.registerExclusivePropGroup(
            e0WherePanel, (e0Where0, e0Where1, e0Where2), e0WherePanel.title(),
            'e0Method')

        layoutWE = qt.QHBoxLayout()
        layoutWE.setContentsMargins(0, 0, 0, 0)
        labelE0 = qt.QLabel('E\u2080 =')
        layoutWE.addWidget(labelE0)
        valueE0 = qt.QLabel('')
        self.registerStatusLabel(valueE0, 'e0', textFormat='.3f')
        layoutWE.addWidget(valueE0)
        layoutWE.addStretch()
        layoutE.addLayout(layoutWE)

        calibPanel = qt.QGroupBox(self)
        calibPanel.setTitle('need energy calibration')
        calibPanel.setCheckable(True)
        self.registerPropWidget(calibPanel, calibPanel.title(),
                                'needECalibration', transformNames='make chi')
        layoutCP = qt.QVBoxLayout()
        layoutCP.setContentsMargins(10, 2, 2, 2)
        layoutCP.setSpacing(2)

        layoutEref0 = qt.QHBoxLayout()
        layoutEref0.setContentsMargins(0, 0, 0, 0)
        howEref0 = qt.QRadioButton('assign to E\u2080 reference E', calibPanel)
        layoutEref0.addWidget(howEref0)
        self.comboEref = qt.QComboBox()
        self.comboEref.setEditable(True)
        energies = XAFSmass.read_energies()
        self.comboEref.addItems(energies)
        self.comboEref.lineEdit().setText('')
        layoutEref0.addWidget(self.comboEref)
        layoutCP.addLayout(layoutEref0)
        self.registerPropWidget(self.comboEref, 'reference energy', 'eRef',
                                convertType=self.energy_selected)

        layoutEref1 = qt.QHBoxLayout()
        layoutEref1.setContentsMargins(0, 0, 0, 0)
        howEref1 = qt.QRadioButton('shift E\u2080 by', calibPanel)
        layoutEref1.addWidget(howEref1)
        eShiftBox = qt.QDoubleSpinBox()
        eShiftBox.setMinimum(-1e5)
        eShiftBox.setMaximum(1e5)
        eShiftBox.setSingleStep(0.001)
        eShiftBox.setDecimals(3)
        eShiftBox.setAccelerated(True)
        self.registerPropWidget(eShiftBox, 'E0 shift', 'eShift')
        self.registerStatusLabel(eShiftBox, 'eShift')
        layoutEref1.addWidget(eShiftBox)
        layoutEref1.addStretch()
        layoutCP.addLayout(layoutEref1)

        layoutEref2 = qt.QHBoxLayout()
        layoutEref2.setContentsMargins(0, 0, 0, 0)
        labelApply = qt.QLabel('apply calibration as ')
        layoutEref2.addWidget(labelApply)
        eShiftKinds = qt.QComboBox()
        eShiftKinds.addItems(xtr.MakeChi.eShiftKinds)
        self.registerPropWidget(eShiftKinds, 'eShift kind', 'eShiftKind')
        layoutEref2.addWidget(eShiftKinds)
        layoutEref2.addStretch()
        layoutCP.addLayout(layoutEref2)

        calibPanel.setLayout(layoutCP)
        layoutE.addWidget(calibPanel)
        self.registerExclusivePropGroup(
            calibPanel, (howEref0, howEref1), calibPanel.title(),
            'eCalibrationMethod')

        e0Panel.setLayout(layoutE)
        layout.addWidget(e0Panel)

        preedgePanel = qt.QGroupBox(self)
        preedgePanel.setFlat(False)
        preedgePanel.setStyleSheet(
            'QGroupBox[flat="false"] {font-weight: bold;}')
        preedgePanel.setTitle('pre-edge background')
        layoutP = gco.QVBoxLayoutAbove()
        layoutP.setContentsMargins(10, 2, 2, 2)

        checkBoxShowPreEdge = qt.QCheckBox('show subtracted', preedgePanel)
        checkBoxShowPreEdge.setChecked(self.properties['subtract_preedge'])
        checkBoxShowPreEdge.toggled.connect(self.subtractPreedgeSlot)
        layoutP.addExtraWidget(checkBoxShowPreEdge)

        self.preedgeRangeWidget = RangeWidgetPre(
            self, plot, 'energy range',
            'fractionMin, fractionMax of range [eMin, E0]', 'pre-edge',
            "#008b8b", "{0[0]:.2f}, {0[1]:.2f}",
            xtr.MakeChi.defaultParams['preedgeWhere'])
        self.registerPropWidget(self.preedgeRangeWidget, 'pre-edge range',
                                'preedgeWhere')
        layoutP.addWidget(self.preedgeRangeWidget)

        self.preedgeStateButtons = gco.StateButtons(
            self, 'exponents', (-4, -3, 0, 1), default=0)
        self.registerPropWidget(self.preedgeStateButtons, 'pre-edge exponents',
                                'preedgeExps')
        layoutP.addWidget(self.preedgeStateButtons)

        preedgePanel.setLayout(layoutP)
        layout.addWidget(preedgePanel)

        postedgePanel = qt.QGroupBox(self)
        postedgePanel.setFlat(False)
        postedgePanel.setTitle('post-edge bckgnd for normalization')
        postedgePanel.setStyleSheet(
            'QGroupBox[flat="false"] {font-weight: bold;}')
        layoutPo = gco.QVBoxLayoutAbove()
        layoutPo.setContentsMargins(10, 2, 2, 2)
        layoutPo.setSpacing(2)

        checkBoxShowPost = qt.QCheckBox('plot it', postedgePanel)
        checkBoxShowPost.setChecked(self.properties['show_post'])
        checkBoxShowPost.toggled.connect(partial(self.showSlot, 'show_post'))
        layoutPo.addExtraWidget(checkBoxShowPost)

        self.postedgeRangeWidget = RangeWidgetPost(
            self, plot, 'energy range', '[min, max] +E0 (eV)',
            'post-edge', "#8b8b00", "{0[0]:.1f}, {0[1]:.1f}",
            xtr.MakeChi.defaultParams['postedgeWhere'])
        self.registerPropWidget(self.postedgeRangeWidget, 'post-edge range',
                                'postedgeWhere')
        layoutPo.addWidget(self.postedgeRangeWidget)

        self.postedgeStateButtons = gco.StateButtons(
            self, 'exponents', (-2, -1, 0, 1, 2), default=0)
        self.registerPropWidget(
            self.postedgeStateButtons, 'post-edge exponents', 'postedgeExps')
        layoutPo.addWidget(self.postedgeStateButtons)

        layoutSt = qt.QHBoxLayout()
        layoutSt.setContentsMargins(2, 2, 2, 2)
        labelSt = qt.QLabel('edge height =')
        layoutSt.addWidget(labelSt)
        valueSt = qt.QLabel('')
        self.registerStatusLabel(valueSt, 'edgeJump', textFormat='.4f')
        layoutSt.addWidget(valueSt)
        layoutSt.addStretch()
        layoutPo.addLayout(layoutSt)

        layoutN = qt.QHBoxLayout()
        layoutN.setContentsMargins(0, 0, 0, 0)
        if not self.properties['subtract_preedge']:
            self.properties['normalize'] = False
        self.checkBoxNormalize = qt.QCheckBox('show edge height normalized')
        self.checkBoxNormalize.setChecked(self.properties['normalize'])
        self.checkBoxNormalize.toggled.connect(self.normalizeSlot)
        layoutN.addWidget(self.checkBoxNormalize)
        if not self.properties['normalize']:
            self.properties['flatten'] = False
        self.checkBoxFlatten = qt.QCheckBox('show flat')
        self.checkBoxFlatten.setEnabled(self.properties['normalize'])
        self.checkBoxFlatten.setChecked(self.properties['flatten'])
        self.checkBoxFlatten.toggled.connect(partial(self.showSlot, 'flatten'))
        layoutN.addWidget(self.checkBoxFlatten)
        layoutPo.addLayout(layoutN)

        postedgePanel.setLayout(layoutPo)
        layout.addWidget(postedgePanel)

        mu0Panel = qt.QGroupBox(self)
        mu0Panel.setFlat(False)
        # stylesheet also sets bold the inner QGroupBox, solve it with font
        # mu0Panel.setStyleSheet(
        #     'QGroupBox[flat="false"] {font-weight: bold;}')
        font = mu0Panel.font()
        font.setBold(True)
        mu0Panel.setFont(font)
        mu0Panel.setTitle('atomic background µ\u2080')
        layoutM = gco.QVBoxLayoutAbove()
        layoutM.setContentsMargins(10, 2, 2, 2)

        checkBoxShowMu0 = qt.QCheckBox('plot µ\u2080', mu0Panel)
        checkBoxShowMu0.setChecked(self.properties['show_mu0'])
        checkBoxShowMu0.toggled.connect(partial(self.showSlot, 'show_mu0'))
        layoutM.addExtraWidget(checkBoxShowMu0)

        mu0PriorPanel = qt.QGroupBox(self)
        mu0PriorPanel.setFlat(False)
        # stylesheet doesn't work for this inner QGroupBox, solve it with font
        font = mu0PriorPanel.font()
        font.setBold(False)
        mu0PriorPanel.setFont(font)
        # mu0PriorPanel.setStyleSheet('QGroupBox {font-weight: light}')
        mu0PriorPanel.setTitle('µ\u2080 prior')
        layoutMP = gco.QVBoxLayoutAbove()
        layoutMP.setContentsMargins(10, 2, 2, 2)
        layoutMP.setSpacing(2)

        checkBoxShowMu0P = qt.QCheckBox('plot µ\u2080 prior', mu0PriorPanel)
        checkBoxShowMu0P.setChecked(self.properties['show_mu0_prior'])
        checkBoxShowMu0P.toggled.connect(
            partial(self.showSlot, 'show_mu0_prior'))
        layoutMP.addExtraWidget(checkBoxShowMu0P)

        layoutMPV = qt.QHBoxLayout()
        layoutMPV.setContentsMargins(0, 0, 0, 0)
        whiteLine = qt.QCheckBox('white line')
        self.registerPropWidget(
            whiteLine, whiteLine.text(), 'mu0PriorIncludeWhiteLine')
        layoutMPV.addWidget(whiteLine)
        layoutMPV.addStretch()
        labelScale = qt.QLabel('×')
        layoutMPV.addWidget(labelScale)
        scaleV = qt.QDoubleSpinBox()
        scaleV.setMinimum(0.4)
        scaleV.setMaximum(2)
        scaleV.setSingleStep(0.01)
        scaleV.setDecimals(2)
        scaleV.setAccelerated(True)
        self.registerPropWidget(scaleV, 'mu0 prior scaling', 'mu0PriorVScale')
        layoutMPV.addWidget(scaleV)
        layoutMPV.addStretch()
        labelSmooth = qt.QLabel('smoothing')
        layoutMPV.addWidget(labelSmooth)
        smoothV = qt.QSpinBox()
        smoothV.setMinimum(0)
        self.registerPropWidget(
            smoothV, 'mu0 prior smoothing', 'mu0PriorSmoothN')
        layoutMPV.addWidget(smoothV)
        layoutMP.addLayout(layoutMPV)

        mu0PriorPanel.setLayout(layoutMP)
        layoutM.addWidget(mu0PriorPanel)

        layoutMPM = qt.QHBoxLayout()
        layoutMPM.setContentsMargins(0, 0, 0, 0)
        labelMethod = qt.QLabel('method')
        layoutMPM.addWidget(labelMethod)
        methods = qt.QComboBox()
        methods.addItems(xtr.MakeChi.mu0methods)
        self.registerPropWidget(methods, 'mu0 method', 'mu0method')
        layoutMPM.addWidget(methods)
        layoutMPM.addStretch()
        layoutM.addLayout(layoutMPM)

        mu0page0 = qt.QWidget()
        layoutMu0page0 = qt.QVBoxLayout()
        layoutMu0page0.setContentsMargins(0, 0, 0, 0)
        mu0page0.setLayout(layoutMu0page0)

        layoutMK = qt.QHBoxLayout()
        nKnots = qt.QLabel('knots')
        self.knotsBox = qt.QSpinBox()
        self.knotsBox.setMinimum(3)
        self.knotsBox.setMaximum(30)
        self.registerPropWidget(self.knotsBox, 'knots, w', 'mu0knots')
        layoutMK.addWidget(nKnots)
        layoutMK.addWidget(self.knotsBox)

        labelMu0kw = qt.QLabel('exp weight')
        self.mu0kpowBox = qt.QSpinBox()
        self.mu0kpowBox.setMinimum(0)
        self.mu0kpowBox.setMaximum(5)
        self.registerPropWidget(self.mu0kpowBox, 'mu0kpow', 'mu0kpow')
        layoutMK.addStretch()
        layoutMK.addWidget(labelMu0kw)
        layoutMK.addWidget(self.mu0kpowBox)
        layoutMK.addStretch()
        layoutMu0page0.addLayout(layoutMK)
        layoutMu0page0.addStretch()

        mu0page1 = qt.QWidget()
        layoutMu0page1 = qt.QVBoxLayout()
        layoutMu0page1.setContentsMargins(0, 0, 0, 0)
        mu0page1.setLayout(layoutMu0page1)

        layoutMSM = qt.QHBoxLayout()
        layoutMSM.setContentsMargins(0, 0, 0, 0)
        labelSm = qt.QLabel('smoothing factor')
        layoutMSM.addWidget(labelSm)
        sm = gco.EDoubleSpinBox(strFormat='{0:.1e}')
        sm.setValue(100)
        sm.setMinimum(1)
        sm.setMaximum(1e10)
        sm.setAccelerated(True)
        self.registerPropWidget(
            sm, 'mu0 smoothing factor', 'mu0smoothingFactor')
        layoutMSM.addWidget(sm)
        layoutMSM.addStretch()

        layoutMu0page1.addLayout(layoutMSM)
        layoutMu0page1.addStretch()

        self.stackedMu0MethodWidget = qt.QStackedWidget()
        self.stackedMu0MethodWidget.addWidget(mu0page0)
        self.stackedMu0MethodWidget.addWidget(mu0page1)
        methods.activated.connect(self.stackedMu0MethodWidget.setCurrentIndex)
        self.stackedMu0MethodWidget.setSizePolicy(
            qt.QSizePolicy.Expanding, qt.QSizePolicy.Fixed)

        layoutM.addWidget(self.stackedMu0MethodWidget)
        mu0Panel.setLayout(layoutM)
        layout.addWidget(mu0Panel)

        layout.addStretch()
        self.setLayout(layout)
        self.wasNeverPlotted = True

    def subtractPreedgeSlot(self, value):
        if len(csi.selectedItems) == 0:
            return
        self.checkBoxNormalize.setEnabled(value)
        self.properties['subtract_preedge'] = value
        if not value:
            self.properties['normalize'] = False
            self.checkBoxNormalize.setChecked(False)

        plot = self.node.widget.plot
        elim = plot.getXAxis().getLimits()
        ylim = plot.getYAxis().getLimits()
        eMiddle = sum(elim) * 0.5
        data = csi.selectedItems[0]
        try:
            pe = interp1d(data.e, data.pre_edge, assume_sorted=True)
            peMiddle = pe(eMiddle)
            ylim = [yl-peMiddle if value else yl+peMiddle for yl in ylim]
            plot.getYAxis().setLimits(*ylim)
        except (ValueError, TypeError):
            pass
        csi.model.needReplot.emit(False, True, 'subtractPreedgeSlot')

    def energy_selected(self, txt=None):
        if txt is None:
            txt = self.comboEref.lineEdit().text()
        calibEnergy = txt
        try:
            calibEnergy = float(txt)
        except ValueError:
            st = str(txt).strip().split()
            if len(st) == 4:
                e = float(st[3])
                e = int(np.floor(e)) if e == np.floor(e) else e
                calibEnergy = e
                self.comboEref.lineEdit().setText(str(calibEnergy))
        return calibEnergy

    def normalizeSlot(self, value):
        if len(csi.selectedItems) == 0:
            return
        self.properties['normalize'] = value
        self.checkBoxFlatten.setEnabled(value)
        plot = self.node.widget.plot
        ylim = list(plot.getYAxis().getLimits())
        data = csi.selectedItems[0]
        ylim[1] *= 1./data.edge_step if value else data.edge_step
        plot.getYAxis().setLimits(*ylim)
        csi.model.needReplot.emit(False, True, 'normalizeSlot')

    def showSlot(self, prop, value):
        self.properties[prop] = value
        csi.model.needReplot.emit(False, True, 'showSlot')

    def extraSetUIFromData(self):
        if len(csi.selectedItems) == 0:
            return
        data = csi.selectedItems[0]
        dtparams = data.transformParams
        self.stackedMu0MethodWidget.setCurrentIndex(dtparams['mu0method'])

        self.checkBoxDerivativeRef.setChecked(dtparams['useERefCurve'])
        if hasattr(data, 'eref'):
            self.checkBoxDerivativeRef.setEnabled(data.eref is not None)

    def extraPlot(self):
        plot = self.node.widget.plot
        for data in csi.allLoadedItems:
            if not self.node.widget.shouldPlotItem(data):
                for extraLine in self.extraLines:
                    legend = '{0}.{1}'.format(data.alias, extraLine)
                    plot.remove(legend, kind='curve')
                continue
            dtparams = data.transformParams
            z = 1 if data in csi.selectedItems else 0

            legend = '{0}.pre_edge'.format(data.alias)
            if not self.properties['subtract_preedge']:
                curve = plot.getCurve(legend)
                if curve is None:
                    plot.addCurve(
                        data.e, data.pre_edge, **self.preEdgePlotParams,
                        color=data.color, z=z, legend=legend, resetzoom=False)
                else:
                    curve.setData(data.e, data.pre_edge)
                    curve.setZValue(z)
            else:
                plot.remove(legend, kind='curve')

            legend = '{0}.post_edge'.format(data.alias)
            if self.properties['show_post']:
                pe = data.post_edge - data.pre_edge if \
                    self.properties['subtract_preedge'] else \
                    np.array(data.post_edge)
                if self.properties['normalize']:
                    if self.properties['flatten']:
                        pe /= data.post_edge - data.pre_edge
                    else:
                        pe /= data.edge_step
                curve = plot.getCurve(legend)
                if curve is None:
                    plot.addCurve(
                        data.e, pe, **self.postEdgePlotParams,
                        color=data.color, z=z, legend=legend, resetzoom=False)
                else:
                    curve.setData(data.e, pe)
                    curve.setZValue(z)
            else:
                plot.remove(legend, kind='curve')

            legend = '{0}.mu0'.format(data.alias)
            if self.properties['show_mu0']:
                pe = data.mu0 - data.pre_edge if \
                    self.properties['subtract_preedge'] else \
                    np.array(data.mu0)
                if self.properties['normalize']:
                    if self.properties['flatten']:
                        pe /= data.post_edge - data.pre_edge
                    else:
                        pe /= data.edge_step
                curve = plot.getCurve(legend)
                if curve is None:
                    plot.addCurve(
                        data.e, pe, **self.mu0PlotParams,
                        color=data.color, z=z, legend=legend, resetzoom=False)
                else:
                    curve.setData(data.e, pe)
                    curve.setZValue(z)
            else:
                plot.remove(legend, kind='curve')

            legend = '{0}.mu0prior'.format(data.alias)
            if self.properties['show_mu0_prior']:
                pe = data.mu0prior - data.pre_edge if \
                    self.properties['subtract_preedge'] else \
                    np.array(data.mu0prior)
                if self.properties['normalize']:
                    if self.properties['flatten']:
                        pe /= data.post_edge - data.pre_edge
                    else:
                        pe /= data.edge_step
                curve = plot.getCurve(legend)
                if curve is None:
                    plot.addCurve(
                        data.e, pe, **self.mu0priorPlotParams,
                        color=data.color, z=z, legend=legend, resetzoom=False)
                else:
                    curve.setData(data.e, pe)
                    curve.setZValue(z)
            else:
                plot.remove(legend, kind='curve')

            legend = '{0}.e0'.format(data.alias)
            if self.properties['show_E0'] and hasattr(data, 'mu'):
                de = data.mu - data.pre_edge if \
                    self.properties['subtract_preedge'] else data.mu
                if self.properties['normalize']:
                    de = np.array(de) / data.edge_step  # make a copy
                # ylim = plot.getYAxis().getLimits()
                ylim = (0, de.max()) if self.properties['subtract_preedge'] \
                    else (de.min(), de.max())
                dy = ylim[1] - ylim[0]
                ylim = ylim[0] + 0.05*dy, ylim[1] - 0.05*dy
                curve = plot.getCurve(legend)
                if curve is None:
                    plot.addCurve(
                        [data.e0, data.e0], ylim, **self.e0PlotParams,
                        color=data.color, z=z, legend=legend, resetzoom=False)
                else:
                    curve.setData([data.e0, data.e0], ylim)
                    curve.setZValue(z)
            else:
                plot.remove(legend, kind='curve')

            legend = '{0}.mu_der'.format(data.alias)
            if self.properties['show_derivative']:
                if dtparams['useERefCurve'] and data.eref_der is not None:
                    mu_der = data.eref_der / data.edge_step \
                        if self.properties['normalize'] else data.eref_der
                else:
                    mu_der = data.mu_der / data.edge_step \
                        if self.properties['normalize'] else data.mu_der
                if mu_der is not None:
                    curve = plot.getCurve(legend)
                    if curve is None:
                        plot.addCurve(data.e, mu_der, **self.derPlotParams,
                                      color=data.color, z=z, legend=legend,
                                      resetzoom=self.wasNeverPlotted)
                    else:
                        curve.setData(data.e, mu_der)
                        curve.setZValue(z)
                    self.wasNeverPlotted = False
                    plot.setGraphYLabel(label="$µd'$", axis='right')
            else:
                plot.remove(legend, kind='curve')

    def extraPlotTransform(self, dataItem, xName, x, yName, y):
        if yName == 'mu':
            if self.properties['subtract_preedge']:
                if self.properties['normalize']:
                    if self.properties['flatten']:
                        return x, dataItem.flat
                    else:
                        return x, dataItem.norm
                else:
                    return x, y-dataItem.pre_edge
            else:
                return x, y
        else:
            return x, y

    @classmethod
    def cursorPositionCallback(cls, label, x, y):
        if len(csi.selectedItems) == 0:
            return '---'
        data = csi.selectedItems[0]

        if label == 'E':
            return x
        elif label == 'k':
            if not (hasattr(data, 'e') and hasattr(data, 'e0')):
                return '---'
            return ((x - data.e0)*xtr.eV2revA)**0.5 if x > data.e0 else '---'
        else:
            return y


class ChiWidget(PropWidget):
    u"""
    Help page under construction

    .. image:: _images/mickey-rtfm.gif
        :width: 309

    test link: `MAX IV Laboratory <https://www.maxiv.lu.se/>`_

    """

    name = u'rebin and make χ'

    properties = {'show_zero_grid_line': True, 'k_range_visible': True,
                  'show_ft_window': True}

    captions = 'pre-edge', 'edge', 'post-edge', 'EXAFS'
    deltas = (['dE', 1.0, 0.1, 10, 0.1],  # label, value, min, max, step
              ['dE', 0.2, 0.02, 2, 0.02],
              ['dE', 0.4, 0.04, 4, 0.02],
              ['dk', 0.025, 0.005, 0.2, 0.001])
    splitters = (['E0+', -20, -200, 0, 1],  # label, value, min, max, step
                 ['E0+', 50., 0, 200, 2],
                 ['kmin', 2.1, 0, 5, 0.1],
                 ['kmax', 'inf', 0, 'inf', 0.1])
    defaultRegions = captions, deltas, splitters
    ftWindowPlotParams = {'linewidth': 0.75, 'linestyle': '-',
                          'color': '#00000044'}

    cursorLabels = ['k', 'E', u'χ']

    def __init__(self, parent=None, node=None):
        super().__init__(parent, node)
        plot = self.node.widget.plot

        layout = qt.QVBoxLayout()

        rebinPanel = qt.QGroupBox(self)
        rebinPanel.setFlat(False)
        rebinPanel.setTitle('rebin energy axis')
        rebinPanel.setCheckable(True)
        rebinPanel.setStyleSheet(
            'QGroupBox[flat="false"] {font-weight: bold;}')
        self.registerPropWidget(rebinPanel, rebinPanel.title(), 'rebinNeeded')
        layoutRB = qt.QVBoxLayout()
        layoutRB.setContentsMargins(2, 2, 2, 2)
        self.regionsWidget = DataRebinWidget(self, self.defaultRegions)
        layoutRB.addWidget(self.regionsWidget)
        # self.acceptButton = qt.QPushButton('Accept regions')
        # layoutRB.addWidget(self.acceptButton, 1)
        self.regionsWidget.regionsChanged.connect(self.acceptRebinRegions)
        self.registerPropWidget(
            [self.regionsWidget.splittersView, self.regionsWidget.deltasView],
            'rebinRegions', 'rebinRegions')
        rebinPanel.setLayout(layoutRB)
        layout.addWidget(rebinPanel)

        krangePanel = qt.QGroupBox(self)
        krangePanel.setFlat(False)
        krangePanel.setTitle('k range')
        krangePanel.setStyleSheet(
            'QGroupBox[flat="false"] {font-weight: bold;}')

        layoutK = qt.QVBoxLayout()
        layoutK.setContentsMargins(2, 0, 2, 2)
        self.krange = RangeWidgetSplit(
            self, plot, ('k min', 'k max'), [0, 10, 0.1, 3],
            [0, None, 0.1, 3], 'k-range', '#ff5500', [1.5, None],
            addVisibleCB=True)
        self.registerPropWidget(self.krange, 'k-range', 'krange')
        layoutK.addWidget(self.krange)
        self.krange.setRangeVisible(self.properties['k_range_visible'])

        layoutDK = self.krange.layout()
        labelDataKMax = qt.QLabel('data k max')
        layoutDK.addWidget(labelDataKMax, 1, 3)
        datakmaxValue = qt.QLabel('')
        self.registerStatusLabel(datakmaxValue, 'datakmax', textFormat='.3f')
        layoutDK.addWidget(datakmaxValue, 1, 4)
        datakmaxButton = qt.QToolButton()
        datakmaxButton.setFixedSize(16, 16)
        datakmaxButton.setIcon(icons.getQIcon('last'))
        datakmaxButton.setToolTip("Set k max from data")
        datakmaxButton.clicked.connect(self.setkmaxFromData)
        color = qt.QColor(gco.COLOR_LOAD_CAN)
        color.setAlphaF(0.32)
        datakmaxButton.setStyleSheet("QToolButton{border-radius: 8px;}"
                                     "QToolButton:hover{background-color: " +
                                     color.name(qt.QColor.HexArgb) + ";}")
        layoutDK.addWidget(datakmaxButton, 1, 5)
        layoutDK.setColumnStretch(6, 1)

        dk = qt.QLabel('dk')
        layoutDK.addWidget(dk, 3, 0, qt.Qt.AlignRight)
        self.dkBox = qt.QDoubleSpinBox()
        self.dkBox.setMinimum(0.005)
        self.dkBox.setMaximum(0.1)
        self.dkBox.setSingleStep(0.005)
        self.dkBox.setDecimals(3)
        self.registerPropWidget(self.dkBox, 'dk', 'dk', dataItems="all")
        layoutDK.addWidget(self.dkBox, 3, 1)

        krangePanel.setLayout(layoutK)
        layout.addWidget(krangePanel)

        layoutkw = qt.QHBoxLayout()
        labelkw = qt.QLabel('×k<sup>w</sup>, w=')
        layoutkw.addWidget(labelkw)
        self.kw = qt.QSpinBox()
        self.kw.setToolTip(u'×k^w, 0 ≤ w ≤ 3')
        self.kw.setMinimum(0)
        self.kw.setMaximum(3)
        self.registerPropWidget((self.kw, labelkw), 'w in k^w', 'kw')
        layoutkw.addWidget(self.kw)
        labelkw2 = qt.QLabel('(applied to all)')
        layoutkw.addWidget(labelkw2)
        layoutkw.addStretch()
        layout.addLayout(layoutkw)

        ftPanel = qt.QGroupBox(self)
        ftPanel.setFlat(False)
        ftPanel.setTitle('FT window')
        ftPanel.setStyleSheet('QGroupBox[flat="false"] {font-weight: bold;}')

        # layoutFT = qt.QVBoxLayout()
        layoutFT = gco.QVBoxLayoutAbove()
        layoutFT.setContentsMargins(2, 0, 2, 2)

        checkBoxShowWindow = qt.QCheckBox('plot it', ftPanel)
        checkBoxShowWindow.setChecked(self.properties['show_ft_window'])
        checkBoxShowWindow.toggled.connect(
            partial(self.showSlot, 'show_ft_window'))
        layoutFT.addExtraWidget(checkBoxShowWindow)

        self.ftWindowKind = qt.QComboBox()
        self.ftWindowKind.addItems(uft.ft_windows)
        self.registerPropWidget(self.ftWindowKind, 'FT window', 'ftWindowKind',
                                dataItems="all", indexToValue=uft.ft_windows,
                                # transformNames='make FT'
                                )
        self.ftWindowKind.currentIndexChanged.connect(self.updateFTwindow)
        layoutFT.addWidget(self.ftWindowKind)

        self.ftWidthAndMin = RangeWidgetFTWidthAndMin(
            self, plot, ('width', 'minimum'), [0., 10.0, 0.1, 2],
            [0.0, 1.0, 0.01, 3], 'FT window\nwidth,FT window\nminimum',
            self.ftWindowPlotParams['color'], [0.5, 0.1])
        self.registerPropWidget(
            self.ftWidthAndMin, 'FT window width and min', 'ftWindowProp',
            dataItems="all",
            # transformNames='make FT'
            )
        layoutFT.addWidget(self.ftWidthAndMin)

        ftPanel.setLayout(layoutFT)
        layout.addWidget(ftPanel)

        layout.addStretch()
        self.setLayout(layout)

    def acceptRebinRegions(self):
        regions = self.regionsWidget.getRegions()
        self.updateProp('rebinRegions', regions)

    def updateStatusWidgets(self):
        if len(csi.selectedItems) == 0:
            return
        super().updateStatusWidgets()
        data = csi.selectedItems[0]
        dtparams = data.transformParams
        self.regionsWidget.setBinNumbers(0, dtparams['nbinOriginal'])
        self.regionsWidget.setBinNumbers(1, dtparams['nbinNew'])
        self.regionsWidget.setBinNumbers(2, dtparams['binDistrNew'])

    def extraSetUIFromData(self):
        if len(csi.selectedItems) == 0:
            return
        data = csi.selectedItems[0]
        dtparams = data.transformParams
        self.regionsWidget.setRegions(dict(dtparams['rebinRegions']))

    def extraPlotActionAfterTransform(self, props):
        if 'kw' in props:
            plot = self.node.widget.plot
            plot.removeCurve('FT window')
            plot.resetZoom()
            nextNodeInd = list(csi.nodes.keys()).index(self.node.name) + 1
            nextNodeName = list(csi.nodes.keys())[nextNodeInd]
            nextNode = csi.nodes[nextNodeName]
            w = nextNode.widget.transformWidget
            plot = nextNode.widget.plot
            plot.resetZoom()
            w.showNegativeSlot(w.properties['show_negative'])

    def setkmaxFromData(self):
        if len(csi.selectedItems) == 0:
            return
        data = csi.selectedItems[0]
        dtparams = data.transformParams
        kmin = dtparams['krange'][0]
        kmax = dtparams['datakmax']
        self.krange.setRange([kmin, kmax])
        self.krange.roi.setMax(kmax)
        self.updateProp('krange', [kmin, kmax])

    def extraPlot(self):
        if len(csi.selectedItems) == 0:
            return
        data = csi.selectedItems[0]
        dtparams = data.transformParams

        plot = self.node.widget.plot
        self.ftWidthAndMin.fromSpinBox(100)
        if hasattr(data, 'ftwindow'):
            legend = 'hline'
            xlim = plot.getXAxis().getLimits()
            if self.properties['show_zero_grid_line']:
                plot.addCurve(xlim, [0, 0], color='gray', legend=legend,
                              resetzoom=False)

            legend = 'FT window'
            if self.properties['show_ft_window']:
                ymax = plot.getYAxis().getLimits()[1] * \
                    RangeWidgetFTWidthAndMin.plotFactor
                plot.addCurve(
                    data.k, data.ftwindow*ymax, yaxis='left',
                    **self.ftWindowPlotParams, legend=legend, resetzoom=False)
            else:
                plot.remove(legend, kind='curve')

    #     for data in csi.allLoadedItems:
    #         if not self.node.widget.shouldPlotItem(data):
    #             continue
    #         z = 1 if data in csi.selectedItems else 0
    #         legend = '{0}.bft'.format(data.alias)
    #         plot.addCurve(
    #             data.k, data.bft,
    #             color='red', z=z, legend=legend, resetzoom=False)
        kw = dtparams['kw']
        if kw == 1:
            skw = u'·k (Å'+u"\u207B"+u'¹)'
        elif kw == 2:
            skw = u'·k² (Å'+u"\u207B"+u'²)'
        elif kw == 3:
            skw = u'·k³ (Å'+u"\u207B"+u'³)'
        else:
            skw = u''
        plot.plotLeftYLabel = u'χ'+skw
        plot.setGraphYLabel(label=plot.plotLeftYLabel, axis='left')
        # plot.setGraphYLabel(label='FT window', axis='right')  # not good

    def updateFTwindow(self, ind):
        self.ftWidthAndMin.setVisible(ind > 1)
        if self.ftWidthAndMin.roi is not None:
            self.ftWidthAndMin.roi.setVisible(ind > 1)

    def showSlot(self, prop, value):
        self.properties[prop] = value
        self.ftWidthAndMin.roi.setVisible(value)
        csi.model.needReplot.emit(False, True, 'showSlot')

    @classmethod
    def cursorPositionCallback(cls, label, x, y):
        if len(csi.selectedItems) == 0:
            return '---'
        data = csi.selectedItems[0]

        if label == 'k':
            return x
        elif label == 'E':
            if not (hasattr(data, 'e') and hasattr(data, 'e0')):
                return '---'
            return data.e0 + x**2/xtr.eV2revA if x > 0 else '---'
        else:
            return y


class FTWidget(PropWidget):
    u"""
    Help page under construction

    .. image:: _images/mickey-rtfm.gif
        :width: 309

    test link: `MAX IV Laboratory <https://www.maxiv.lu.se/>`_

    """

    name = 'make FT'

    properties = {'show_negative': False, 'show_Re': False, 'show_Im': True}

    extraLines = '-|ft|', 'ft.re', 'ft.im'
    ftrPlotParams = {'linewidth': 0.7, 'linestyle': '-.'}
    ftiPlotParams = {'linewidth': 0.7, 'linestyle': ':'}

    def __init__(self, parent=None, node=None):
        super().__init__(parent, node)
        layout = qt.QVBoxLayout()

        layoutMR = qt.QHBoxLayout()
        layoutMR.setContentsMargins(0, 0, 0, 0)
        labelRMax = qt.QLabel('r max')
        layoutMR.addWidget(labelRMax)
        rmax = qt.QDoubleSpinBox()
        rmax.setMinimum(0)
        rmax.setMaximum(40)
        rmax.setSingleStep(0.01)
        rmax.setDecimals(2)
        rmax.setAccelerated(True)
        self.registerPropWidget(rmax, 'r max', 'rmax')
        layoutMR.addWidget(rmax)
        layoutMR.addStretch()
        layout.addLayout(layoutMR)

        forceFT0 = qt.QCheckBox('force FT(0)=0')
        self.registerPropWidget(forceFT0, forceFT0.text(), 'forceFT0')
        layout.addWidget(forceFT0)

        self.checkBoxShowNegative = qt.QCheckBox('show negative part')
        self.checkBoxShowNegative.setChecked(self.properties['show_negative'])
        self.checkBoxShowNegative.toggled.connect(self.showNegativeSlot)
        layout.addWidget(self.checkBoxShowNegative)

        self.checkBoxShowRe = qt.QCheckBox('show FT.Re')
        self.checkBoxShowRe.setChecked(self.properties['show_Re'])
        self.checkBoxShowRe.toggled.connect(partial(self.showSlot, 'show_Re'))
        layout.addWidget(self.checkBoxShowRe)

        self.checkBoxShowIm = qt.QCheckBox('show FT.Im')
        self.checkBoxShowIm.setChecked(self.properties['show_Im'])
        self.checkBoxShowIm.toggled.connect(partial(self.showSlot, 'show_Im'))
        layout.addWidget(self.checkBoxShowIm)

        layout.addStretch()
        self.setLayout(layout)

    def showNegativeSlot(self, value):
        self.properties['show_negative'] = value
        plot = self.node.widget.plot
        _, ymax = plot.getYAxis().getLimits()
        ylim = [-ymax, ymax] if value else [0, ymax]
        plot.getYAxis().setLimits(*ylim)
        csi.model.needReplot.emit(False, True, 'showNegativeSlot')

    def showSlot(self, prop, value):
        self.properties[prop] = value
        csi.model.needReplot.emit(False, True, 'showSlot')

    def extraPlot(self):
        if len(csi.selectedItems) == 0:
            return
        plot = self.node.widget.plot
        for data in csi.allLoadedItems:
            if not self.node.widget.shouldPlotItem(data):
                for extraLine in self.extraLines:
                    legend = '{0}.{1}'.format(data.alias, extraLine)
                    plot.remove(legend, kind='curve')
                continue
            z = 1 if data in csi.selectedItems else 0

            legend = '{0}.-|ft|'.format(data.alias)
            if self.properties['show_negative']:
                plotProps = dict(data.plotProps[self.node.name]['ft'])
                curve = plot.getCurve(legend)
                if curve is None:
                    plot.addCurve(
                        data.r, -data.ft, color=data.color, z=z, legend=legend,
                        resetzoom=False, **plotProps)
                else:
                    curve.setData(data.r, -data.ft)
                    curve.setZValue(z)
            else:
                plot.remove(legend, kind='curve')

            legend = '{0}.ft.re'.format(data.alias)
            if self.properties['show_Re']:
                if self.properties['show_negative']:
                    y = data.ftr
                else:
                    y = np.array(data.ftr)
                    y[y < 0] = 0
                curve = plot.getCurve(legend)
                if curve is None:
                    plot.addCurve(
                        data.r, y, **self.ftrPlotParams,
                        color=data.color, z=z, legend=legend, resetzoom=False)
                else:
                    curve.setData(data.r, y)
                    curve.setZValue(z)
            else:
                plot.remove(legend, kind='curve')

            legend = '{0}.ft.im'.format(data.alias)
            if self.properties['show_Im']:
                if self.properties['show_negative']:
                    y = data.fti
                else:
                    y = np.array(data.fti)
                    y[y < 0] = 0
                curve = plot.getCurve(legend)
                if curve is None:
                    plot.addCurve(
                        data.r, y, **self.ftiPlotParams,
                        color=data.color, z=z, legend=legend, resetzoom=False)
                else:
                    curve.setData(data.r, y)
                    curve.setZValue(z)
            else:
                plot.remove(legend, kind='curve')

        data = csi.selectedItems[0]
        dtparams = data.transformParams
        kw = dtparams['kw']
        if kw == 1:
            skw = u'·k'
            skf = u'²'
        elif kw == 2:
            skw = u'·k²'
            skf = u'³'
        elif kw == 3:
            skw = u'·k³'
            skf = u'⁴'
        else:
            skw = u''
            skf = u'¹'
        plot.plotLeftYLabel = u'FT[χ'+skw+u'] (Å'+u"\u207B"+skf+u')'
        plot.setGraphYLabel(label=plot.plotLeftYLabel, axis='left')
