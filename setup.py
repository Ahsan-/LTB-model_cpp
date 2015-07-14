from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [
    Extension("emcee_wrapper",#extra_compile_args = ['-O3'], 
              sources=["cmb.pyx"],
              language="c++",
              )
   ]
)

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [
    Extension("scipy_tools",#extra_compile_args = ['-O3'], 
              sources=["scipy_tools.pyx"],
              language="c++",
              )
   ]
)
