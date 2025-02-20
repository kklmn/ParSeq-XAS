# -*- coding: utf-8 -*-
u"""
Data nodes
----------

Experimental signals
~~~~~~~~~~~~~~~~~~~~

One design goal of ParSeq is the ability of data import at *any* data node, not
only at the head node(s) of a pipeline. Nonetheless, a few nodes of the
ParSeq-XAS pipeline are considered to be *input* nodes -- those listed in this
seqtion -- where experimental signals are loaded. The nodes :class:`NodeIF` and
:class:`NodeIE` do the same job and can recieve single or multiple fluorescence
or electron yield signals. The idea of the splitting into two separate nodes is
to reduce the need for data format re-definition when switching between
multi-channel and single-channel signals. The node :class:`NodeIXES` receives
a 2D intensity array of a HERFD scan.

.. note::
    Data format definitions are explained in .... In particular, any data
    channel can be defined by a Python expression. Therefore, if 'eref' does
    not exist as an array in an hdf5 file or a column in a column file, it is
    possible to define it as e.g. `np.log(Col5/Col4)`.

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
    The following three arrays are required:
      | 'eraw': the original energy axis in eV; it will later transform to 'e',
      | 'i0': the I₀ signal -- instensity upstream of the sample,
      | 'itr': the transmitted intensity.

    An optional array 'eref' can be given to be used for energy calibration.
    If used, it should contain an absorption spectrum of a foil.
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
    The following three arrays are required:
      | 'eraw': the original energy axis in eV; it will later transform to 'e',
      | 'i0': the I₀ signal -- instensity upstream of the sample,
      | 'ify': the partial fluorescence yield signal.

    An optional array 'eref' can be given to be used for energy calibration.
    If used, it should contain an absorption spectrum of a foil.

    .. note::
        Data format definitions are explained in .... In particular,
        multiple fluorescence channels can either be loaded as a sum or a
        collection of individual spectra. The data format can be defined by a
        Python list or a list comprehension expression, e.g
        `[f'Col{n}' for n in range(10, 17)]`. One can also insert it into the
        `sum()` function.
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
    The following three arrays are required:
      | 'eraw': the original energy axis in eV; it will later transform to 'e',
      | 'i0': the I₀ signal -- instensity upstream of the sample,
      | 'iey': the partial fluorescence yield signal.

    An optional array 'eref' can be given to be used for energy calibration.
    If used, it should contain an absorption spectrum of a foil.
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
    The following three arrays are required:
      | 'eraw': the original energy axis in eV; it will later transform to 'e',
      | 'i0': the I₀ signal -- instensity upstream of the sample,
      | 'xes2D': the 2D intensity array in (scan axis (DCM energy) vs
                 meridional detector pixel) coordinates.

    An optional array 'eref' can be given to be used for energy calibration.
    If used, it should contain an absorption spectrum of a foil.
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
    The following two arrays are required:
      | 'eraw': the original energy axis in eV; it will transform to 'e',
      | 'muraw': the original absortion coefficient; it will transform to 'mu'.

    An optional array 'eref' can be given to be used for energy calibration.
    If used, it should contain an absorption spectrum of a foil.
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
    The following two arrays are required:
      | 'k': the photoelectron wavenumber axis as an equidistant mesh,
      | 'chi': the EXAFS function multiplied by kʷ.
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
    The following four arrays are required:
      | 'r': the uncorrected distance axis as an equidistant mesh,
      | 'ft': the absolute Fourier-transformed EXAFS function multiplied by kʷ.
      | 'ftr': the real Fourier-transformed EXAFS function multiplied by kʷ.
      | 'fti': the imag Fourier-transformed EXAFS function multiplied by kʷ.
    """

    name = 'FT, χ(r)'
    description = "EXAFS function in r-space"
    icon = "doc/_images/icon-xas-ft.png"
    ref = "nodes.html#fourier-transformed-exafs-function-r"
    arrays = OrderedDict()
    arrays['r'] = dict(qUnit=u'Å', role='x', plotLabel=r'$r$')
    arrays['ft'] = dict(
        qLabel=u'|FT|', qUnit=u'Å\u207B\u00B9', role='yleft',
        plotLabel=r'|FT($\chi$)|', plotUnit=r'Å$^{-1}$',
        plotParams=dict(linewidth=1.5, linestyle='-'))
    arrays['ftr'] = dict(
        qLabel=u'Re(FT)', qUnit=u'Å\u207B\u00B9', role='yleft',
        plotLabel=r'Re(FT($\chi$))', plotUnit=r'Å$^{-1}$',
        plotParams=dict(linewidth=0.7, linestyle='-.', hidden=True))
    arrays['fti'] = dict(
        qLabel=u'Im(FT)', qUnit=u'Å\u207B\u00B9', role='yleft',
        plotLabel=r'Im(FT($\chi$))', plotUnit=r'Å$^{-1}$',
        plotParams=dict(linewidth=0.7, linestyle=':', hidden=True))


class NodeBFT(cno.Node):
    """
    The following two arrays are required:
      | 'bftk': the photoelectron wavenumber axis as an equidistant mesh,
      | 'bft': the EXAFS function multiplied by kʷ.
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
