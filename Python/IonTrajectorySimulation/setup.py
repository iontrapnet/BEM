#from distutils.core import setup
#from distutils.extension import Extension
from setuptools import setup, Extension
from Cython.Distutils import build_ext

ext_modules=[
    Extension("simulation",
              ["simulation.pyx"],
              #libraries=["m"] # Unix-like specific
              )
]
setup(
  name = 'Simulation',
  cmdclass = {"build_ext": build_ext},
  ext_modules = ext_modules
)