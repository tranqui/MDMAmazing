#!/usr/bin/env python3

from setuptools import setup, Extension, find_packages
from glob import glob
import os

cpp_args = ['-std=c++14', '-mtune=native', '-march=native', '-Ofast', '-Wall', '-Wextra']
link_args =  ['-Wl,--unresolved-symbols=report-all']

ext_source = []
for root,dirs,files in os.walk('src/pybind11'):
    ext_source += ['%s/%s' % (root,f) for f in files]
ext_path = [path.replace('src/pybind11', 'mdma').replace('.cc','') for path in ext_source]

ext_modules = [
    Extension(module_path, [source_path],
              include_dirs=['include', 'pybind11/include', '/usr/include/eigen3'],
              language='c++',
              extra_compile_args = cpp_args,
              extra_link_args = link_args)
    for module_path, source_path in zip(ext_path, ext_source)
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
    package_dir={'': 'src'},
    packages=find_packages('src'),
    ext_modules=ext_modules,
    install_requires=['sphinx', 'pybind11', 'numpy', 'scipy', 'pandas', 'mpi4py', 'natsort', 'progressbar', 'lxml', 'beautifulsoup4']
 )
