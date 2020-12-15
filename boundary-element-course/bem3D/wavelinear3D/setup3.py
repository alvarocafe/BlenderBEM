# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 21:53:35 2017

@author: Admin
"""

from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name = "integracao3",
  cmdclass = {"build_ext": build_ext},
  ext_modules =
  [
    Extension("integracao3",
              ["integracao3.pyx"],
              extra_compile_args = ["-O3"],
              )
  ]
)