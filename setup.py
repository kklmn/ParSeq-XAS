# -*- coding: utf-8 -*-
from setuptools import setup

import os.path as osp
import codecs

__dir__ = osp.abspath(osp.dirname(__file__))


def read(pathnames):
    with codecs.open(osp.join(__dir__, *pathnames), 'r') as fp:
        return fp.read()


def get_version():
    inLines = read(('parseq_XAS', 'version.py')).splitlines()
    for line in inLines:
        if line.startswith('__versioninfo__'):
            versioninfo = eval(line[line.find('=')+1:])
            version = '.'.join(map(str, versioninfo))
            return version
    else:
        raise RuntimeError("Unable to find version string.")


long_description = u"""
EXAFS analysis pipeline
=======================

A pipeline for the [ParSeq framework](https://github.com/kklmn/ParSeq) that
implements data processing of EXAFS spectra.

This pipeline also serves as an example for creating analysis nodes,
transformations that connect these nodes and widgets that set options and
parameters of the transformations.

Dependencies
------------

* [ParSeq](https://github.com/kklmn/ParSeq) -- the framework package,
* [silx](https://github.com/silx-kit/silx) -- used for plotting and Qt imports.

How to use
----------

Either install ParSeq and this pipeline application by their installers or put
their folders `parseq` and `parseq_XAS` near by (i.e. in the same folder) and
run `python XAS_start.py --help` to see the accepted options. Load a ready
project from `saved` folder from the GUI or from the starting command line.
"""

setup(
    name='parseq_XAS',
    version=get_version(),
    description='A pipeline for data processing of XAS spectra',
    long_description=long_description,
    long_description_content_type='text/x-rst',
    author='Konstantin Klementiev',
    author_email='konstantin.klementiev@gmail.com',
    project_urls={'Source': 'https://github.com/kklmn/ParSeq-XAS'},
    platforms='OS Independent',
    license='MIT License',
    keywords='data-analysis pipeline framework gui synchrotron spectroscopy',
    # python_requires=,
    zip_safe=False,  # True: build zipped egg, False: unzipped
    packages=['parseq_XAS'],
    package_data={
        'parseq_XAS': ['data/*.*', 'data/foils/*.*', 'data/HERFD/*.*',
                       'doc/_images/*.*', 'saved/*.*']},
    scripts=['parseq_XAS/XAS_start.py'],
    install_requires=['numpy>=1.8.0', 'scipy>=1.10.0', 'matplotlib>=2.0.0',
                      'sphinx>=1.6.2', 'h5py', 'silx>=1.1.0', 'hdf5plugin'],
    classifiers=['Development Status :: 5 - Production/Stable',
                 'Intended Audience :: Science/Research',
                 'Natural Language :: English',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'License :: OSI Approved :: MIT License',
                 'Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering',
                 'Topic :: Software Development',
                 'Topic :: Software Development :: User Interfaces']
    )
