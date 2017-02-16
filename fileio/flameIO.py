
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
from flame import Machine
from .fileIO import LatticeFile

class FlameLatticeFile(LatticeFile):
    element_type=['SOURCE','DRIFT','SBEND','QUADRUPOLE','SOLENOID','RFCAVITY','STRIPPER','EDIPOLE','EQUAD', 'BPM' ,'MARKER','GENERIC']
    def __init__(self, filename='', read_from=''):
        LatticeFile.__init__(self)
        self.filename=filename
        if read_from != '':
            self.parseFrom(read_from)

    def parseFrom(self, lattice_file_name):
        with open(lattice_file_name, 'r') as f:
            self.machine=Machine(f)
        lattice_dict=self.machine.conf()
        self.IonEk = lattice_dict['IonEk']
        self.IonEs = lattice_dict['IonEs']
        self.IonChargeStates = lattice_dict['IonChargeStates']
        self.lattice = lattice_dict['elements']
        self.HdipoleFitMode=float(lattice_dict.get("HdipoleFitMode",1.0))
        self.MpoleLevel=float(lattice_dict.get("MpoleLevel",2.0))
        self.NCharge=lattice_dict['NCharge']

        n_charge_state=len(self.IonChargeStates)
        self.BaryCenter=[]
        self.SMatrix=[]
        for i in range(n_charge_state):
            temp='BaryCenter{}'.format(i)
            self.BaryCenter.append(lattice_dict[temp])
            temp = 'S{}'.format(i)
            self.SMatrix.append(lattice_dict[temp])


        for ele in self.lattice:
            self.addElement(ele['name'],ele['type'], **ele)
            self.appendToBeamline('MAIN_LINE', ele['name'].upper())

    def checkType(self, typename, parameterName=None):
        #Should be checked by the flame parser
        return True



