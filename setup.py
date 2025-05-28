import sys

from pybind11 import get_cmake_dir
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages

__version__ = "0.0.1"

ext_modules = [
    Pybind11Extension("bisr_dpw_cpp_routines",
        ["cpp/cpp_routines.cpp"],
        define_macros = [('VERSION_INFO', __version__)],
        ),
]

setup(
    name="rna_barrier_bisr",
    version=__version__,
    install_requires=["networkx",],
    packages=["rna_barrier_bisr"],
    cmdclass={"build_ext":build_ext},
    ext_modules=ext_modules
)
