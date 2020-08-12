#!/usr/bin/env python3

from setuptools import setup, Extension, find_packages
from glob import glob
import sys, os

required_modules = ['numpy', 'scipy', 'pandas', 'natsort', 'progressbar', 'lxml', 'beautifulsoup4']
if '--with-mpi' in sys.argv:
    required_modules += ['mpi4py']

requirements = []
for module in required_modules:
    try: exec("import %s" % module)
    except: requirements += [module]

ext_modules = []
if '--with-pybind11' in sys.argv:
    import pybind11
    sys.argv.remove('--with-pybind11')
    requirements += ['pybind11']

    cpp_args = ['-std=c++14', '-mtune=native', '-march=native', '-Ofast', '-Wall', '-Wextra']
    link_args =  []

    ext_source = []
    for root,dirs,files in os.walk('src/pybind11'):
        ext_source += ['%s/%s' % (root,f) for f in files]
        ext_path = [path.replace('src/pybind11', 'mdma').replace('.cc','') for path in ext_source]

    ext_modules = [
        Extension(module_path, [source_path],
                  include_dirs=['include', pybind11.get_include(False), pybind11.get_include(True), '/usr/include/eigen3'],
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
    version='1.0.1a',
    package_dir={'': 'src'},
    packages=find_packages('src'),
    ext_modules=ext_modules,
    install_requires=requirements
 )
