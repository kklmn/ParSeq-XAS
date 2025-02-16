# -*- coding: utf-8 -*-
u"""
Data nodes
----------

Experimental signals
~~~~~~~~~~~~~~~~~~~~

.. autoclass:: NodeIT

.. autoclass:: NodeIF

.. autoclass:: NodeIE

.. autoclass:: NodeIXES

Absorption coefficient µ(E)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: NodeMu

EXAFS function χ(k)
~~~~~~~~~~~~~~~~~~~

.. autoclass:: NodeChi

Fourier-transformed EXAFS function χ(r)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: NodeFT

Fourier-filtered EXAFS function χ\u0303(k)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: NodeBFT

"""

__author__ = "Konstantin Klementiev"
__date__ = "15 Feb 2025"

import sys; sys.path.append('..')  # analysis:ignore
from collections import OrderedDict
# import hdf5plugin  # needed to prevent h5py's "OSError: Can't read data"

from parseq.core import nodes as cno


class NodeIT(cno.Node):
    """
    Work in progress
    """
    name = 'currents Tr'
    description = "transmission signals"
    icon = "doc/_images/icon-xas-cur-tr.png"
    ref = "nodes.html#experimental-signals"
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
    """
    Work in progress
    """
    name = 'counts PFY'
    description = "partial fluorescence signals"
    icon = "doc/_images/icon-xas-cur-FY.png"
    ref = "nodes.html#experimental-signals"
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
    """
    Work in progress
    """
    name = 'currents TEY/TFY'
    description = "total electron yield or total fluorescence signals"
    icon = "doc/_images/icon-xas-cur-TFY.png"
    ref = "nodes.html#experimental-signals"
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
    """
    Work in progress
    """
    name = '2D XES'
    description = "HERFD maps: meridional pixel as x and incident energy as y"
    icon = "doc/_images/icon-xas-cur-HERFD.png"
    ref = "nodes.html#experimental-signals"
    arrays = OrderedDict()
    arrays['eraw'] = dict(qLabel='E', qUnit='eV', role='y', plotLabel=r'$E$')
    arrays['i0'] = dict(  # not plotted, therefore role='1D'
        qLabel='I0', qUnit='A', role='1D')
    arrays['xes2D'] = dict(
        qLabel='XES2D', qUnit='counts', role='2D',
        plotLabel=['tangential pixel', 'eraw'])
    arrays['eref'] = dict(role='optional', qLabel='Eref')


class NodeMu(cno.Node):
    """
    Work in progress
    """
    name = u'µd'
    description = "µd (optical thickness)"
    icon = "doc/_images/icon-xas-mu.png"
    ref = "nodes.html#absorption-coefficient-e"
    arrays = OrderedDict()
    arrays['e'] = dict(qLabel='E', qUnit='eV', role='x', plotLabel=r'$E$',
                       raw='eraw')
    arrays['mu'] = dict(
        qLabel=u'µd', role='yleft', plotLabel=r'$\mu d$', raw='muraw',
        plotParams=dict(linewidth=1.5, linestyle='-'))
    arrays['eref'] = dict(role='optional', qLabel='Eref')


class NodeChi(cno.Node):
    """
    Work in progress
    """
    name = u'χ(k)'
    description = "EXAFS function in k-space"
    icon = "doc/_images/icon-xas-chi.png"
    ref = "nodes.html#exafs-function-k"
    arrays = OrderedDict()
    arrays['k'] = dict(
        qUnit=u'Å\u207B\u00B9', role='x', plotLabel=r'$k$',
        plotUnit=r'Å$^{-1}$')
    arrays['chi'] = dict(
        qLabel=u'χ', role='yleft', plotLabel=r'$\chi$',
        plotParams=dict(linewidth=1.5, linestyle='-'))


class NodeFT(cno.Node):
    """
    Work in progress
    """
    name = 'FT, χ(r)'
    description = "EXAFS function in r-space"
    icon = "doc/_images/icon-xas-ft.png"
    ref = "nodes.html#fourier-transformed-exafs-function-r"
    arrays = OrderedDict()
    arrays['r'] = dict(qUnit=u'Å', role='x', plotLabel=r'$r$')
    arrays['ft'] = dict(
        qLabel=u'|FT(χ)|', qUnit=u'Å\u207B\u00B9', role='yleft',
        plotLabel=r'|FT($\chi$)|', plotUnit=r'Å$^{-1}$',
        plotParams=dict(linewidth=1.5, linestyle='-'))


class NodeBFT(cno.Node):
    """
    Work in progress
    """
    name = 'BFT, χ\u0303(k)'
    description = "FT filtered EXAFS function in k-space"
    icon = "doc/_images/icon-xas-bft.png"
    ref = "nodes.html#fourier-filtered-exafs-function-k"
    arrays = OrderedDict()
    arrays['bftk'] = dict(
        qUnit=u'Å\u207B\u00B9', role='x', plotLabel=r'$k$',
        plotUnit=r'Å$^{-1}$')
    arrays['bft'] = dict(
        qLabel=u'χ\u0303\u00A0', role='yleft',
        plotParams=dict(linewidth=1.5, linestyle='-'))
