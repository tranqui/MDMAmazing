Getting started
###############

Prerequisites
=============

python modules
---------------

Prerequisites:

* numpy
* scipy

Install via pip::

  pip install numpy scipy

as superuser to install these system-wide or::

  pip install --user numpy scipy

to install without elevated privileges.

Installation
============

The repository is hosted at `<https://github.com/tranqui/MDMAmazing>`_.

Clone the repository and install it via::

  git clone https://github.com/tranqui/MDMAmazing>
  cd MDMAmazing
  ./setup.py build
  ./setup.py --user install

Optional: Installing LAMMPS
===========================

To use the wrapper for LAMMPS, it must be compiled with as a shared library and PyLammps must be installed.
Follow the instructions  `here <https://lammps.sandia.gov/doc/Howto_pylammps.html#system-wide-installation>`_ for details of this process.
