#!/usr/bin/env python3

class PotentialBase:
    def __init__(self, natoms, temperature, composition=None, d=3, **kwargs):
        """If composition is None then it is either trivial (single component) or
        must be specified with an initial snapshot."""
        self.natoms = natoms
        self.temperature = temperature
        self.composition = composition
        self.d = d

        if 'density' in kwargs and 'pressure' in kwargs:
            raise ValueError('must specify (at most) one of density or pressure (NVT vs NPT simulations)!')
        elif 'density' in kwargs: self.density = kwargs['density']
        elif 'pressure' in kwargs: self.pressure = kwargs['pressure']

    @property
    def units(self):
        raise NotImplementedError("must specify simulation units!")

    @property
    def num_species(self):
        raise NotImplementedError("must specify number of atomic species!")

    def initialise_potential(self, simulation):
        raise NotImplementedError("must initialise potential!")

    @property
    def timestep(self):
        raise NotImplementedError("must specify a timestep!")

    @property
    def initial_thermalisation_steps(self):
        raise NotImplementedError("must specify timestep for initial thermalisation!")

    @property
    def thermostat_damping_interval(self):
        raise NotImplementedError("must specify damping inverval for thermostat!")

class KobAnderson(PotentialBase):
    def __init__(self, *args, **kwargs):
        if 'composition' not in kwargs:
            kwargs['composition'] = [0.8, 0.2]
        super().__init__(*args, **kwargs)

    @property
    def units(self):
        return 'lj'

    @property
    def num_species(self):
        return 2

    def initialise_potential(self, simulation):
        for i in range(self.num_species): simulation.mass(i+1, 1.0)
        simulation.pair_style('lj/cut', 2.5)
        simulation.pair_coeff(1, 1, 1.0, 1.00, 2.5)
        simulation.pair_coeff(1, 2, 1.5, 0.80, 2.5)
        simulation.pair_coeff(2, 2, 0.5, 0.88, 2.5)

    @property
    def timestep(self):
        return 0.005

    @property
    def initial_thermalisation_steps(self):
        return 1000

    @property
    def thermostat_damping_interval(self):
        return 5
