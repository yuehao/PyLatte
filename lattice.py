from . import element
from .fileio import fileIO

import numpy as np
import copy


class Lattice(object):
    def __init__(self):
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
        self.dimension = 6
        self.radiation = True
        self.useline=lattice_file.useline
        if self.useline=='':
            self.useline=lattice_file.beamlineList[-1]['NAME']

        self.make_sequence()


    def make_sequence(self):
        pass



    def floor_map(self):
        pass




