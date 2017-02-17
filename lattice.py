from . import element
from .fileio import fileIO

import numpy as np
import copy


class Lattice(object):
    _modes=['tpsa', 'tracking']
    _radiation=['ISR','eloss']
    _particle=['electron', 'proton', 'position', 'antiproton', 'ion']

    def __init__(self, lattice_file, mode='TPSA', ref_energy=1e9, particle='electron', dimension=6, radiation='ISR', charge=-1, mass=510998.9):
        if isinstance(lattice_file, fileIO.LatticeFile):
            if len(lattice_file.beamlineList) > 0:
                self.lattice_file = lattice_file
            else:
                print("No beamline defined, exiting\n")
                exit(-1)
        else:
            print("No lattice information given, exiting\n")
            exit(-1)

        self.ref_energy = ref_energy
        self.dimension = dimension
        self.radiation = radiation
        self.charge = charge
        self.mass = mass

        self.useline=lattice_file.useline
        if self.useline=='':
            self.useline=lattice_file.beamlineList[-1]['NAME']

        self.make_sequence()


    def make_sequence(self):
        pass



    def floor_map(self):
        pass




