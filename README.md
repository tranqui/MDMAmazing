# MDMAmazing
Molecular Dynamics Management and Analysis


---

## Prerequisites

Eigen3, pybind11 and openmpi must be installed.

On ubuntu install Eigen3 and openmpi via

    sudo apt install libeigen3-dev
    sudo apt-get install libopenmpi-dev

Then follow the instructions on https://pybind11.readthedocs.io/en/stable/basics.html to install pybind11.

## Installation

From inside the code repository execute

    python setup.py install

as admin to install system-wide. Alternatively, if you do not have elevated rights privileges you can install for the user via

    python setup.py install --user
