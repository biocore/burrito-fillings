#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, biocore development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

__version__ = '0.1.1'

from setuptools import find_packages, setup
from distutils.command.build_py import build_py

classes = """
    Development Status :: 1 - Planning
    License :: OSI Approved :: BSD License
    Topic :: Software Development :: Libraries
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python
    Programming Language :: Python :: 2.7
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

long_description = """The burrito-fillings project"""

setup(name='burrito-fillings',
      cmdclass={'build_py': build_py},
      version=__version__,
      license='BSD',
      description=\
        'burrito-fillings: burrito application controllers for bioinformatics',
      long_description=long_description,
      author="biocore",
      author_email="gregcaporaso@gmail.com",
      maintainer="biocore",
      maintainer_email="gregcaporaso@gmail.com",
      url='https://github.com/biocore/burrito-fillings',
      packages=find_packages(),
      install_requires=['scikit-bio >= 0.2.1, < 0.3.0', 'burrito  < 1.0.0'],
      classifiers=classifiers)
