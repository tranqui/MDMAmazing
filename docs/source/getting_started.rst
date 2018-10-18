Getting started
###############

.. warning:: I have only tried installing this on my own (GNU/Linux) machine, so some of these steps may fail on other systems particularly on Windows or Mac OS X operating systems. If you encounter any problems or have any suggestions please contact the `author <index.html#author>`_.

Prerequisites
=============

Python 3
--------

MDMAmazing is written exclusively in Python 3. As such Python 3 must be installed on your system beforehand. Most GNU/Linux systems should come with Python 3 preinstalled, otherwise you can obtain Python from either of:

  * The official source: `<https://www.python.org/downloads/>`_
  * Through `anaconda <https://www.anaconda.com/download>`_ which offers a version of Python prebuilt with libraries for scientific computing.

Installation
============

The repository is hosted at `<https://github.com/tranqui/MDMAmazing>`_.

Clone the repository and install it via::

  git clone https://github.com/tranqui/MDMAmazing>
  cd MDMAmazing
  python3 setup.py build
  python3 setup.py --user install

Optional: Installing LAMMPS
===========================

To use the wrapper for LAMMPS, it must be compiled as a shared library and PyLammps must be installed.
Follow the instructions  `here <https://lammps.sandia.gov/doc/Howto_pylammps.html#system-wide-installation>`_ for details of this process.
