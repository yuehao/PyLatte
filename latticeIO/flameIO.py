
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import copy
from flame import Machine
from .fileIO import LatticeFile

class FlameLatticeFile(LatticeFile):
    element_type=['SOURCE','DRIFT','SBEND','QUADRUPOLE','SOLENOID','RFCAVITY','STRIPPER','EDIPOLE','EQUAD', 'BPM' ,'MARKER','GENERIC']
    def __init__(self, read_from=''):
        LatticeFile.__init__(self)
        self.lattice_format='Flame'
        self.machine = None
        self.linename=None
        if read_from != '':
            self.parseFrom(read_from)

    def isDrift(self, ele_name):
        parent_type = self.getParentElements(ele_name)[-1]
        if 'DRIFT' in parent_type:
            temp=self.getElementProperties(ele_name)
            return temp
        else:
            return False

    def isDipole(self, ele_name):
        parent_type = self.getParentElements(ele_name)[-1]
        if 'BEND' in parent_type or 'DIPOLE' in parent_type:
            temp=self.getElementProperties(ele_name)
            temp.update({'TYPE':parent_type})
            if 'K1' not in temp:
                temp['K1']=0
            return temp
        else:
            return False

    def isQuadrupole(self, ele_name):
        parent_type = self.getParentElements(ele_name)[-1]
        if 'QUAD' in parent_type:
            temp=self.getElementProperties(ele_name)
            if 'B2' in temp:
                temp['K1']=temp['B2']
            elif 'V' in temp:
                    temp['K1'] = temp['V']
            else:
                temp['K1']=0
            return temp
        else:
            return False
    def isSolenoid(self, ele_name):
        parent_type = self.getParentElements(ele_name)[-1]
        if 'SOLENOID' in parent_type:
            temp=self.getElementProperties(ele_name)
            return temp
        else:
            return False
    def isCavity(self, ele_name):
        parent_type = self.getParentElements(ele_name)[-1]
        if 'RFCAVITY' in parent_type:
            temp=self.getElementProperties(ele_name)
            return temp
        else:
            return False


    def parseFrom(self, lattice_file_name):
        with open(lattice_file_name, 'rb') as f:
            self.machine=Machine(f)
        lattice_dict=self.machine.conf()
        self.lattice = copy.deepcopy(lattice_dict['elements'])
        self.linename = lattice_dict['name']
        '''['AMU',
         'Eng_Data_Dir',
         'HdipoleFitMode',
         'IonChargeStates',
         'IonEk',
         'IonEs',
         'IonW',
         'IonZ',
         'MpoleLevel',
         'NCharge',

         'Stripper_IonChargeStates',
         'Stripper_NCharge',

         'sim_type']
        '''

        self.AMU = lattice_dict['AMU']
        self.Eng_Data_Dir = lattice_dict['Eng_Data_Dir']
        self.HdipoleFitMode = float(lattice_dict.get("HdipoleFitMode", 1.0))
        self.IonChargeStates = lattice_dict['IonChargeStates']
        self.IonEk = lattice_dict['IonEk']
        self.IonEs = lattice_dict['IonEs']
        self.IonW = lattice_dict.get('IonW', -1)
        self.IonZ = lattice_dict.get('IonZ', -1)
        self.MpoleLevel = float(lattice_dict.get("MpoleLevel", 2.0))
        self.NCharge=lattice_dict['NCharge']
        self.Stripper_IonChargeStates = lattice_dict['Stripper_IonChargeStates']
        self.Stripper_NCharge=lattice_dict['Stripper_NCharge']
        self.sim_type = lattice_dict['sim_type']

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
            self.appendToBeamline(self.linename, ele['name'].upper())

        self.setUseLine(self.linename)

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
        self.result_ref=np.zeros((len(self.flame_result), 10))
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
            self.result_ref[ind, 8] = itm[1].ref_IonQ
            self.result_ref[ind, 9] = itm[1].ref_IonZ


            self.result_moments[ind, 0] = ind_ele
            self.result_moments[ind, 1] = itm[1].pos
            self.result_moments[ind, 2:8] = itm[1].moment0_env[0:6]
            self.result_moments[ind, 8:] = itm[1].moment0_rms[0:6]

    def _output(self, f):
        f.write()


    def write(self, out_file_name):
        with open(out_file_name, 'w') as f:
            self._output(f)

    def plot_references(self, axis, axis_layout, ploting='ENERGY'):
        if isinstance(ploting, str):
            plist=[ploting.upper(), None]
        elif isinstance(ploting, list) and len(ploting)==2:
            plist=[ploting[0].upper(), ploting[1].upper()]

        colors=['b','g']
        for i in range(2):
            if i==0:
                ax=axis
            elif plist[i] is None:
                break
            else:
                axis_m = axis.twinx()
                ax=axis_m
            if plist[i]=='E' or plist[i]=='ENERGY':
                ax.plot(self.result_ref[:, 1], self.result_ref[:, 5]/1e6, color=colors[i],label='Kinetic Energy')
                ax.set_ylabel("Energy (MeV)")


            elif plist[i]=='PH' or plist[i]=='PHASE':
                ax.plot(self.result_ref[:, 1], self.result_ref[:, 6],color=colors[i], label='Ref Phase')
                ax.set_ylabel("Phase (rad)")


            elif plist[i]=='BG' or plist[i]=='MOMENTUM':
                ax.plot(self.result_ref[:, 1], self.result_ref[:, 3], color=colors[i], label='Norm. Momentum')
                ax.set_ylabel(r"$\beta\gamma$")

        axis.legend(loc='best')
        axis_layout.set_xlabel("s [m]")
        self.plotBeamline(axis_layout)

    def plot_transverse(self, axis, axis_layout, plot_orbit='xy', plot_size=None):
        colors = ['b', 'g', 'r', 'y']
        ci=0

        if plot_orbit is not None:
            axis.set_ylabel("Centroid (mm)")
            if 'X' in plot_orbit.upper():
                axis.plot(self.result_moments[:, 1], self.result_moments[:, 2], color=colors[ci], label='X centroid')
                ci+=1
            if 'Y' in plot_orbit.upper():
                axis.plot(self.result_moments[:, 1], self.result_moments[:, 4], color=colors[ci], label='Y centroid')
                ci += 1
            axis_m = axis.twinx()
        else:
            axis_m=axis
        if plot_size is not None:
            axis_m.set_ylabel("rms beamsize (mm)")
            if 'X' in plot_size.upper():
                axis_m.plot(self.result_moments[:, 1], self.result_moments[:, 8], color=colors[ci], label='X rms size')
                ci += 1
            if 'Y' in plot_size.upper():
                axis_m.plot(self.result_moments[:, 1], self.result_moments[:, 10], color=colors[ci], label='Y rms size')
                ci += 1
        axis.legend(loc='best')
        axis_layout.set_xlabel("s [m]")
        self.plotBeamline(axis_layout)
            
            












