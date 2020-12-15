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
  name = "integracao2",
  cmdclass = {"build_ext": build_ext},
  ext_modules =
  [
    Extension("integracao2",
              ["integracao2.pyx"],
              extra_compile_args = ["-O3", "-fopenmp"],
              extra_link_args=['-fopenmp']
              )
  ]
)