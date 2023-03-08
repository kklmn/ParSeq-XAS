# -*- coding: utf-8 -*-
__author__ = "Konstantin Klementiev"
__date__ = "28 Jan 2022"

import os.path as osp
from os.path import dirname as up
import sys; sys.path.append('..')  # analysis:ignore

import parseq.core.singletons as csi

dirname = up(osp.abspath(__file__))


def load_test_data():
    fNames = [[osp.join(dirname, 'data', 'Cu_lnt1.fio'), ['3', 'Col5', 6, '']],
              [osp.join(dirname, 'data', 'Cu_lnt2.fio'), [3, 5, 6, '']],
              [osp.join(dirname, 'data', 'Cu_rt1.fio'), [3, 5, 6, '']],
              [osp.join(dirname, 'data', 'Cu_rt2.fio'), [3, 5, 6, '']],
              [osp.join(dirname, 'data', 'Cu2O_lnt1.fio'), [0, 5, 6, '']],
              [osp.join(dirname, 'data', 'Cu2O_lnt2.fio'), [0, 5, 6, '']],
              [osp.join(dirname, 'data', 'CuO_lnt.fio'), [0, 5, 6]]]
    conversionFactors = [None, 'counts', 'counts', None]

    # h5base = "silx:{0}".format(osp.join(dirname, 'data', 'Cu-flyScans.h5'))
    # h5Names = ['entry737', 'entry738']
    # h5Format = [
    #     'measurement/mono1_energy',
    #     'd["measurement/albaem01_ch1"] + d["measurement/albaem01_ch4"]',
    #     'measurement/albaem01_ch2']
    # h5ConversionFactors = [None, 1e-3, 1e-3]

    rootItem = csi.dataRootItem
    rootItem.kwargs['runDownstream'] = True

    group0 = rootItem.insert_item('metal', colorPolicy='loop1')
    dataFormat = dict(dataSource=fNames[0][1], lastSkipRowContains='Col ',
                      conversionFactors=conversionFactors)
    data = [fn[0] for fn in fNames[:4]]
    group0.insert_data(data, dataFormat=dataFormat)

    # another way to make a group, res=list:
    group1, = rootItem.insert_data('oxides', colorPolicy='loop2')
#    group1, = rootItem.insert_data('oxides', color='green')
    dataFormat = dict(dataSource=fNames[4][1], lastSkipRowContains='Col ',
                      conversionFactors=conversionFactors)
    data = [fn[0] for fn in fNames[4:7]]
    group1.insert_data(data, dataFormat=dataFormat)

    # group2 = rootItem.insert_item(
    #     'metal-flyScan', colorPolicy='gradient', color1='red', color2='blue',
    #     colorAutoUpdate=True)
    # data = ['::/'.join([h5base, e]) for e in h5Names]
    # dataFormat = dict(
    #     dataSource=h5Format, conversionFactors=h5ConversionFactors)
    # group2.insert_data(data, dataFormat=dataFormat)

    # fNames2 = ['Cu-WT_EXAFS_001.dat', 'Cu-WT_EXAFS_002.dat']

    # group3 = rootItem.insert_item('FY', colorPolicy='loop1')
    # dataFormat = dict(dataSource=[1, 'd["3"]+d["4"]', list(range(8, 15))],
    #                   lastSkipRowContains='Col ')
    # data = [osp.join(dirname, 'data', fName) for fName in fNames2]
    # group3.insert_data(data, dataFormat=dataFormat,
    #                    originNodeName='currents FY')

    csi.allLoadedItems[:] = []
    csi.allLoadedItems.extend(csi.dataRootItem.get_items())
