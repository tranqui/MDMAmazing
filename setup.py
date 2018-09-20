#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name='MDMAmazing',
    description='Molecular Dynamics Management and Analysis',
    long_description=open('README.md').read(),
    author='Joshua Robinson',
    author_email='joshua.robinson@bristol.ac.uk',
    url='https://github.com/tranqui/MDMAmazing.git',
    license='GNU General Public License v3.0',
    version='0.1dev',
    packages=find_packages(),
    install_requires=['numpy', 'scipy', 'pandas', 'mpi4py', 'natsort', 'progressbar']
 )
