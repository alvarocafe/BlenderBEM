# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 21:53:35 2017

@author: Admin
"""

from distutils.core import setup
from Cython.Build import cythonize
setup(ext_modules=cythonize('cindex2.pyx'))
