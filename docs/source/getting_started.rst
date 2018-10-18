Getting started
###############

.. warning:: I have only tried installing this on my own (GNU/Linux) machine, so some of these steps may fail on other systems particularly on Windows or Mac OS X operating systems. If you encounter any problems or have any suggestions please contact the `author <index.html#author>`_.

Prerequisites
=============

Coming soon.

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
