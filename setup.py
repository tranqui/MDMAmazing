#!/usr/bin/env python3

from setuptools import setup, Extension, find_packages

cpp_args = ['-std=c++14', '-mtune=native', '-march=native', '-Ofast']

ext_modules = [
    Extension('mdma/spatial/distance', ['src/spatial/distance.cc'],
              include_dirs=['pybind11/include', '/usr/include/eigen3'],
              language='c++',
              extra_compile_args = cpp_args,
    ),
]

setup(
    name='MDMAmazing',
    description='Molecular Dynamics Management and Analysis',
    long_description=open('README.md').read(),
    author='Joshua Robinson',
    author_email='joshua.robinson@bristol.ac.uk',
    url='https://github.com/tranqui/MDMAmazing.git',
    license='GNU General Public License v3.0',
    version='0.1.dev0',
    packages=find_packages(),
    #ext_modules=ext_modules,
    install_requires=['pybind11', 'numpy', 'scipy', 'pandas', 'mpi4py', 'natsort', 'progressbar']
 )
