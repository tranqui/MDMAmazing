#!/usr/bin/env python3

import os, string, tempfile, shutil
import numpy as np
from glob import glob
from natsort import natsorted
from itertools import combinations
from jeli import periodic

from progressbar import ProgressBar, Percentage, Bar, ETA
widgets = [Percentage(), ' ', Bar(), ' ', ETA()]

from copy import deepcopy
from lxml import etree as ElementTree
from lxml.etree import Element
import bz2
import iniparser

class MonodisperseHardSpheres:
    size_ratios = np.array([1.])
    def density(phi): return 6 * phi / np.pi
    def volume_fraction(density): return density / (6 * np.pi)

class PaddyPolydisperseHardSpheres:
    # 5-component system has these diameters
    size_ratios = np.array([0.7986,0.8609,0.8993,0.9376,1.0000])
    # Effective volume of a single particle (mean of the sizes, as system is equimolar)
    #effective_volume = np.pi*np.sum(size_ratios**3)/(6*len(size_ratios))

    def density(phi): return 6 * phi / np.pi
    def volume_fraction(density): return density / (6 * np.pi)

# Box size for a system with the above parameters.
def box_size(phi, n, system=MonodisperseHardSpheres):
    volume = n / system.density(phi)
    return volume**(1./3)

class DynamoTrajectory:
    def create(path, ncells, phi, system=MonodisperseHardSpheres):
        os.makedirs(path)

        config_path = '%s/0.xml' % path
        os.system('dynamod -m0 --density %.12f -C %r -o %s > /dev/null 2>&1' % (system.density(phi), ncells, config_path))

        trajectory = DynamoTrajectory(path)
        trajectory.apply_uniform_dispersity(phi, system.size_ratios)
        trajectory.save(config_path)
        return trajectory

    def growth_create(old_path, new_path, target_phi, num_relaxations=1):
        old_dynamo = DynamoTrajectory(old_path)
        start = os.path.abspath(old_dynamo.latest_snapshot)
        end = os.path.abspath('%s/0.xml' % new_path)
        collisions = num_relaxations*old_dynamo.relaxation_collisions
        growth_rate = 1/collisions

        os.makedirs(new_path)
        work_dir = tempfile.mkdtemp(prefix='dynamo_')
        #command = '(cd %s; dynarun %s --engine=3 --growth-rate %.8f --target-pack-frac %.4f -o %s > /dev/null 2>&1)' % (work_dir, start, growth_rate, target_phi, end)
        command = '(cd %s; dynarun %s --engine=3 --growth-rate %.8f --target-pack-frac %.4f -o %s)' % (work_dir, start, growth_rate, target_phi, end)
        os.system(command)
        shutil.rmtree(work_dir, ignore_errors=True)

        trajectory = DynamoTrajectory(new_path)
        return trajectory

    def __init__(self, path, d=3):
        self.d = d
        self.root = path
        self.parse_config()
        self.parse_xml(self.latest_snapshot)

    def parse_config(self):
        config_file = '%s/settings.ini' % self.root
        self.settings = iniparser.ConfigFile(config_file)

        # If empty initialise
        if not self.settings: 
            self.settings = iniparser.ConfigFile()
            relaxation = dict()
            relaxation['equilibrated'] = False
            self.settings['Relaxation'] = relaxation

    def run(self, collisions):
        start = os.path.abspath(self.latest_snapshot)
        end = os.path.abspath(self.new_snapshot)

        work_dir = tempfile.mkdtemp(prefix='dynamo_')
        os.system('(cd %s; dynarun %s -E -c %d -o %s > /dev/null 2>&1)' % (work_dir, start, collisions, end))
        shutil.rmtree(work_dir, ignore_errors=True)

        self.parse_xml(end)

    def run_trajectory(self, nframes, dump_time):
        start = os.path.abspath(self.latest_snapshot)

        work_dir = tempfile.mkdtemp(prefix='dynamo_')
        os.system('(cd %s; dynarun %s -E -f %f --snapshot %f -o dump.xml > /dev/null 2>&1)' % (work_dir, start, (nframes+0.5)*dump_time, dump_time))

        paths = glob('%s/*.xml.bz2' % work_dir)
        paths = natsorted([p for p in paths if 'output' not in p])
        assert len(paths) == nframes

        coordinates = np.empty((nframes, self.n, self.d))
        for frame,path in enumerate(paths):
            if 'output' in path: continue

            f = bz2.BZ2File(path)
            parser = ElementTree.parse(f)

            particles = parser.findall("//Pt/P")
            coordinates[frame] = np.array([list(map(float, p.attrib.values())) for p in particles]).reshape(-1,self.d)

        shutil.rmtree(work_dir, ignore_errors=True)

        return coordinates

    def determine_relaxation_time(self, collisions_per_frame = 500):
        try: return self.settings['Relaxation']['tau_collisions']
        except: pass

        print('finding equilibration time...')
        print('collisions\tself-overlap')

        start = self.coords.copy()
        frames = 0
        Q = 1.
        while Q > np.exp(-1):
            self.run(collisions_per_frame)
            Q = periodic.self_overlap(start, self.coords, self.box)
            frames += 1
            print('\t%d\t%.2f' % (collisions_per_frame * frames, Q))

        # Estimate equilibration time
        tau_collisions = collisions_per_frame * frames
        print('approximate equilibration time:', tau_collisions, 'collisions')
        self.settings['Relaxation']['tau_collisions'] = tau_collisions
        self.save_settings()

        return tau_collisions

    @property
    def relaxation_collisions(self):
        return self.determine_relaxation_time()

    @property
    def equilibrated(self):
        return self.settings['Relaxation']['equilibrated']
    def equilibrate(self, num_relaxations = 100):
        if self.equilibrated: return

        tau_collisions = self.relaxation_collisions

        print('equilibrating over %d relaxation times...' % num_relaxations)
        progress = ProgressBar(widgets=widgets)
        for _ in progress(range(num_relaxations)):
            self.run(tau_collisions)

        self.settings['Relaxation']['equilibrated'] = True
        self.save_settings()

    def save(self, path=None):
        if path is None: path = self.new_snapshot
        self.xml['tree'].write(path, pretty_print=True)
        self.save_settings()

    def save_settings(self):
        self.settings.write('%s/settings.ini' % self.root)

    def clear_snapshots(self):
        for snap in self.snapshots: os.remove(snap)
        self.save('%s/0.xml' % self.root)

    @property
    def snapshots(self):
        return natsorted(glob('%s/*.xml' % self.root))

    @property
    def latest_snapshot(self):
        return self.snapshots[-1]

    @property
    def new_snapshot(self):
        latest = self.latest_snapshot.split('/')[-1]
        latest = int(latest.split('.')[0])
        return '%s/%d.xml' % (self.root,latest+1)

    @property
    def n(self):
        return len(self.xml['particles'])

    @property
    def coords(self):
        try: return self._coords
        except:
            self._coords = np.array([p.find('P').attrib.values() for p in self.xml['particles']], dtype=np.longdouble)
            return self._coords

    @coords.setter
    def coords(self, x):
        self._coords = x
        for particle in self.xml['particles']:
            i = int(particle.attrib['ID'])
            P = particle.find('P').attrib
            for c,key in enumerate(P): P[key] = str(x[i,c])

    @property
    def velocities(self):
        try: return self._velocities
        except:
            self._velocities = np.array([p.find('V').attrib.values() for p in self.xml['particles']], dtype=np.longdouble)
            return self._velocities

    @property
    def xyzr(self):
        return np.concatenate([self.coords,
                               0.5*self.diameters.reshape(-1,1)],axis=1)

    @property
    def coords_trajectory(self):
        trajectory = []
        for snap in self.snapshots:
            self.parse_xml(snap)
            trajectory += [self.coords]
        return trajectory

    @property
    def xyzr_trajectory(self):
        trajectory = []
        for snap in self.snapshots:
            self.parse_xml(snap)
            trajectory += [self.xyzr]
        return trajectory

    @property
    def effective_diameter(self):
        effective_volume = self.volume_fraction / self.density
        return (6 * effective_volume / np.pi)**(1./3)

    def parse_xml(self, path):
        try: del self._coords
        except: pass

        # Read the xml file
        parser = ElementTree.XMLParser(remove_blank_text=True)
        self.xml = {}
        self.xml['tree'] = ElementTree.parse(path, parser)

        self.xml['root'] = self.xml['tree'].getroot()
        self.xml['particles'] = self.xml['root'].find("ParticleData")
        self.xml['simulation'] = self.xml['root'].find('Simulation')
        self.xml['interactions'] = self.xml['simulation'].find('Interactions')
        self.xml['genus'] = self.xml['simulation'].find('Genus')

        # System size information
        self.N = len(self.xml['particles'].getchildren())
        self.box = self.xml['simulation'].find('SimulationSize')
        self.box = np.array([float(self.box.attrib[dim]) for dim in ['x','y','z'][:self.d]])
        self.volume = np.product(self.box)
        self.density = self.N / self.volume

        # Find the diameters of the particles, assuming additive interactions.
        self.diameters = np.empty(self.N)
        for uij in self.xml['interactions']:
            uij_range = uij.find('IDPairRange')
            if uij_range.attrib['Type'] == 'Pair': continue
            if uij_range.attrib['Type'] == 'All':
                self.diameters[:] = float(uij.attrib['Diameter'])
                break

            uij_range = uij_range.getchildren()[0].attrib
            start = int(uij_range['Start'])
            end = int(uij_range['End'])
            self.diameters[start:end+1] = float(uij.attrib['Diameter'])

        # Compute exact volume fraction.
        intrinsic_volumes = np.pi * self.diameters**3 / 6
        self.volume_fraction = np.sum(intrinsic_volumes) / self.volume

    def apply_uniform_dispersity(self, phi, size_ratios):
        names = list(string.ascii_uppercase)[:len(size_ratios)]

        bin_indices = np.array_split(np.arange(self.N), len(size_ratios))
        bin_indices = [(str(bin[0]), str(bin[-1])) for bin in bin_indices]
        bin_widths = [1+int(bin[1])-int(bin[0]) for bin in bin_indices]

        # Volume occupied by each partial component
        partial_weights = [n*s**3  for n,s in zip(bin_widths,size_ratios)]
        partial_volumes = np.pi / 6 * np.array(partial_weights)
        # Effective one-component intrinsic volume
        effective_volume = np.sum(partial_volumes) / self.N

        # Normalise the size distribution so we achieve the correct volume fraction.
        norm = phi / (effective_volume*self.density)
        sizes = size_ratios * norm**(1./3)
        for (a,b),size, in zip(bin_indices,sizes):
            self.diameters[int(a):int(b)+1] = size
        self.volume_fraction = phi

        for species in self.xml['genus']:
            self.xml['genus'].remove(species)
        for name,diameter,(start,end) in zip(names,sizes,bin_indices):
            species = Element("Species", Mass="1", Name=name, Type="Point")
            species.append(Element("IDRange", Type="Ranged", Start=start, End=end))
            self.xml['genus'].append(species)

        for uij in self.xml['interactions']:
            self.xml['interactions'].remove(uij)

        for name,diameter,(start,end) in zip(names,sizes,bin_indices):
            interaction = Element("Interaction", Type="HardSphere", Diameter=str(diameter), Name='%s%sInteraction' % (name,name))
            id_range = Element("IDPairRange", Type="Single")
            id_range.append(Element("IDRange", Type="Ranged", Start=start, End=end))
            interaction.append(id_range)
            self.xml['interactions'].append(interaction)

        for (nameA,diameterA,(startA,endA)), (nameB,diameterB,(startB,endB)) in combinations(zip(names,sizes,bin_indices), 2):
            interaction = Element("Interaction", Type="HardSphere", Diameter=str(0.5*(diameterA+diameterB)), Name='%s%sInteraction' % (nameA,nameB))
            id_range = Element("IDPairRange", Type="Pair")
            id_range.append(Element("IDRange", Type="Ranged", Start=startA, End=endA))
            id_range.append(Element("IDRange", Type="Ranged", Start=startB, End=endB))
            interaction.append(id_range)
            self.xml['interactions'].append(interaction)

        # Permute the particles to randomise the size assignments
        self.coords = self.coords[np.random.permutation(self.N)]
