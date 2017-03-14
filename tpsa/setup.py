from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize(
           "tpsa.pyx",                 # our Cython source
           sources=["tpsalib.cpp"],  # additional source file(s)
           language="c++",             # generate C++ code
      ))


