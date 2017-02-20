from . import element
from .fileio import fileIO

import numpy as np
import copy


class Lattice(object):
    _mode_type=['TPSA', 'TRACKING']
    _radiation_type=['ISR','ELOSS']


    def __init__(self, lattice_file = None, mode='TPSA', ref_energy=1e9, radiation='ISR'):

        self._lattice_file = None
        if lattice_file is not None:
            self.read_lattice(lattice_file)
        self._mode=mode
        self._ref_energy = ref_energy
        self._radiation = radiation
        self.sequence=None

        self.make_sequence()

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self,value):
        if value.upper() in Lattice._mode_type:
            self._mode = value.upper()
        else:
            raise KeyError

    @property
    def radiation(self):
        return self._radiation

    @mode.setter
    def radiation(self, value):
        if value.upper() in Lattice._radiation_type:
            self._radiation = value.upper()
        else:
            raise KeyError

    def read_lattice(self, lattice_file):
        if isinstance(lattice_file, fileIO.LatticeFile):
            if len(lattice_file.beamlineList) > 0:
                self.lattice_file = lattice_file
            else:
                print("No beamline defined, exiting\n")
                exit(-1)
        else:
            print("No lattice information given, exiting\n")
            exit(-1)

        self.useline = lattice_file.useline
        if self.useline == '':
            self.useline = lattice_file.beamlineList[-1]['NAME']
        self._make_sequence()


    def _make_sequence(self):
        self.sequence=[]




    def floor_map(self):
        pass




