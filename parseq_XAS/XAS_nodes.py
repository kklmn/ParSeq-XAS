# -*- coding: utf-8 -*-
__author__ = "Konstantin Klementiev"
__date__ = "13 Feb 2025"

import sys; sys.path.append('..')  # analysis:ignore
from collections import OrderedDict
# import hdf5plugin  # needed to prevent h5py's "OSError: Can't read data"

from parseq.core import nodes as cno


class NodeIT(cno.Node):
    name = 'currents Tr'
    description = "transmission signals"
    icon = "doc/_images/icon-xas-cur-tr.png"
    arrays = OrderedDict()
    arrays['eraw'] = dict(qLabel='E', qUnit='eV', role='x', plotLabel=r'$E$')
    arrays['i0'] = dict(
        qLabel='I0', qUnit='A', role='yleft', plotLabel=r'$I_0$',
        plotParams=dict(linewidth=1.0, linestyle='-'))
    arrays['itr'] = dict(
        qLabel='Itr', qUnit='A', role='yright', plotLabel=r'$I_{\rm tr}$',
        plotParams=dict(linewidth=1.8, linestyle='-'))
    arrays['eref'] = dict(role='optional', qLabel='Eref')


class NodeIF(cno.Node):
    name = 'counts PFY'
    description = "partial fluorescence signals"
    icon = "doc/_images/icon-xas-cur-FY.png"
    arrays = OrderedDict()
    arrays['eraw'] = dict(qLabel='E', qUnit='eV', role='x', plotLabel=r'$E$')
    arrays['i0'] = dict(
        qLabel='I0', qUnit='A', role='yleft', plotLabel=r'$I_0$',
        plotParams=dict(linewidth=0.8, linestyle='-'))
    arrays['ify'] = dict(
        qLabel='PFY', qUnit='counts', role='yright', plotLabel=r'$I_{\rm FY}$',
        plotParams=dict(linewidth=1.5, linestyle='-'))
    arrays['eref'] = dict(role='optional', qLabel='Eref')


class NodeIE(cno.Node):
    name = 'currents TEY/TFY'
    description = "total electron yield or total fluorescence signals"
    icon = "doc/_images/icon-xas-cur-TFY.png"
    arrays = OrderedDict()
    arrays['eraw'] = dict(qLabel='E', qUnit='eV', role='x', plotLabel=r'$E$')
    arrays['i0'] = dict(
        qLabel='I0', qUnit='A', role='yleft', plotLabel=r'$I_0$',
        plotParams=dict(linewidth=0.8, linestyle='-'))
    arrays['iey'] = dict(
        qLabel='TEY', qUnit='counts', role='yright',
        plotLabel=r'$I_{\rm TEY}$',
        plotParams=dict(linewidth=1.5, linestyle='-'))
    arrays['eref'] = dict(role='optional', qLabel='Eref')


class NodeIXES(cno.Node):
    name = '2D XES'
    description = "HERFD maps: meridional pixel as x and incident energy as y"
    icon = "doc/_images/icon-xas-cur-HERFD.png"
    arrays = OrderedDict()
    arrays['eraw'] = dict(qLabel='E', qUnit='eV', role='y', plotLabel=r'$E$')
    arrays['i0'] = dict(  # not plotted, therefore role='1D'
        qLabel='I0', qUnit='A', role='1D')
    arrays['xes2D'] = dict(
        qLabel='XES2D', qUnit='counts', role='2D',
        plotLabel=['tangential pixel', 'eraw'])
    arrays['eref'] = dict(role='optional', qLabel='Eref')


class NodeMu(cno.Node):
    name = u'µd'
    description = "µd (optical thickness)"
    icon = "doc/_images/icon-xas-mu.png"
    arrays = OrderedDict()
    arrays['e'] = dict(qLabel='E', qUnit='eV', role='x', plotLabel=r'$E$',
                       raw='eraw')
    arrays['mu'] = dict(
        qLabel=u'µd', role='yleft', plotLabel=r'$\mu d$', raw='muraw',
        plotParams=dict(linewidth=1.5, linestyle='-'))
    arrays['eref'] = dict(role='optional', qLabel='Eref')


class NodeChi(cno.Node):
    name = u'χ(k)'
    description = "EXAFS function in k-space"
    icon = "doc/_images/icon-xas-chi.png"
    arrays = OrderedDict()
    arrays['k'] = dict(
        qUnit=u'Å\u207B\u00B9', role='x', plotLabel=r'$k$',
        plotUnit=r'Å$^{-1}$')
    arrays['chi'] = dict(
        qLabel=u'χ', role='yleft', plotLabel=r'$\chi$',
        plotParams=dict(linewidth=1.5, linestyle='-'))


class NodeFT(cno.Node):
    name = 'FT, χ(r)'
    description = "EXAFS function in r-space"
    icon = "doc/_images/icon-xas-ft.png"
    arrays = OrderedDict()
    arrays['r'] = dict(qUnit=u'Å', role='x', plotLabel=r'$r$')
    arrays['ft'] = dict(
        qLabel=u'|FT(χ)|', qUnit=u'Å\u207B\u00B9', role='yleft',
        plotLabel=r'|FT($\chi$)|', plotUnit=r'Å$^{-1}$',
        plotParams=dict(linewidth=1.5, linestyle='-'))


class NodeBFT(cno.Node):
    name = 'BFT, χ\u0303(k)'
    description = "FT filtered EXAFS function in k-space"
    icon = "doc/_images/icon-xas-bft.png"
    arrays = OrderedDict()
    arrays['bftk'] = dict(
        qUnit=u'Å\u207B\u00B9', role='x', plotLabel=r'$k$',
        plotUnit=r'Å$^{-1}$')
    arrays['bft'] = dict(
        qLabel=u'χ\u0303\u00A0', role='yleft',
        plotParams=dict(linewidth=1.5, linestyle='-'))
