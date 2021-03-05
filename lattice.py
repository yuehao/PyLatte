from . import element
from .latticeIO import fileIO

import numpy as np
import copy


class Lattice(object):
    _mode_type=['TPSA', 'TRACKING']
    _radiation_type=['ISR','ELOSS']


    def __init__(self, lattice_file = None, lattice_name=None,  mode='TPSA', ref_energy=1e9, radiation='ISR'):


        if lattice_file is not None:
            self.get_lattice(lattice_file, lattice_name)
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

    @radiation.setter
    def radiation(self, value):
        if value.upper() in Lattice._radiation_type:
            self._radiation = value.upper()
        else:
            raise KeyError

    def get_lattice(self, lattice_file, lattice_name):
        self.lattice_file=lattice_file
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
        self.element_position, element_name_list = lattice_file.elementPosInUseLine, lattice_file.useLineList


        for elename in element_name_list:
            propertylist=lattice_file.getElementProperties(elename)











