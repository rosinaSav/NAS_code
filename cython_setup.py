'''
Author: Rosina Savisaar.
Have Cython convert cython_func.pyx into C.
'''

from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'cython-func',
  ext_modules = cythonize("cython_func.pyx"),
)
