
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
        self.machine = None
        if read_from != '':
            self.parseFrom(read_from)

    def isDrift(self, ele, parent_type=None):
        if ele['TYPE']=='DRIFT' or parent_type=='DRIFT':
            return {'L':ele['L']}
        else:
            return None

    def isDipole(self, ele, parent_type=None):
        if ele['TYPE'] == 'SBEND' or parent_type=='SBEND':
            return {'L': ele['L'], 'ANGLE':ele['PHI'], 'K1':0}
        else:
            return None

    def isQuadrupole(self, ele, parent_type=None):
        if ele['TYPE'] == 'QUADRUPOLE' or parent_type=='QUADRUPOLE':
            return {'L': ele['L'], 'K1':ele['B2']}
        else:
            return None
    def isSolenoid(self,ele, parent_type=None):
        if ele['TYPE'] == 'SOLENOID' or parent_type=='SOLENOID':
            return {'L': ele['L'], 'B':ele['B']}
        else:
            return None
    def isCavity(self,ele, parent_type=None):
        if ele['TYPE'] == 'RFCAVITY' or parent_type=='RFCAVITY':
            return {'L': ele['L']}
        else:
            return None


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
            self.setUseLine('MAIN_LINE')

    def checkType(self, typename, parameterName=None):
        #Should be checked by the flame parser
        return True

    def run(self):
        if self.machine is None:
            print('No lattice is defined')
            return
        s=self.machine.allocState({})
        self.machine.propagate(s,0,1)
        startn=1
        stopn=len(self.machine)-1
        self.flame_result=self.machine.propagate(s,startn,stopn-startn+1,observe=range(len(self.machine)))
        self._collect_data()


    def _collect_data(self):
        arg_kws=['pos', 'ref_beta','ref_bg', 'ref_gamma', 'ref_IonEk', 'ref_IonEs',
                 'ref_IonQ', 'ref_IonW', 'ref_IonZ', 'ref_phis', 'ref_SampleIonK',
                 'beta', 'bg', 'gamma', 'IonEk', 'IonEs','IonQ', 'IonW','IonZ', 'phis',
                 'SampleIonK', 'moment0','moment0_rms', 'moment0_env', 'moment1']
        self.result_ref=np.zeros((len(self.flame_result), 8))
        self.result_moments = np.zeros((len(self.flame_result), 14))
        for ind in range(len(self.flame_result)):
            itm=self.flame_result[ind]
            ind_ele=itm[0]
            
            self.result_ref[ind, 0] = ind_ele
            self.result_ref[ind, 1] = itm[1].pos
            self.result_ref[ind, 2] = itm[1].ref_beta
            self.result_ref[ind, 3] = itm[1].ref_bg
            self.result_ref[ind, 4] = itm[1].ref_gamma
            self.result_ref[ind, 5] = itm[1].ref_IonEk
            self.result_ref[ind, 6] = itm[1].ref_phis
            self.result_ref[ind, 7] = itm[1].ref_SampleIonK

            self.result_moments[ind, 0] = ind_ele
            self.result_moments[ind, 1] = itm[1].pos
            self.result_moments[ind, 2:8] = itm[1].moment0_env[0:6]
            self.result_moments[ind, 8:] = itm[1].moment0_rms[0:6]

            
            
            
            
            
            
            












