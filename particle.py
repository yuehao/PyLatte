import numpy as np
import scipy.constants as const



class Particle(object):
    def __init__(self, species):
        self.charge=1
        self.species=species.upper()
        self.mass_eV=0
        self.energy=0
        self.momentum=0
        self._process()
    def _process(self):
        if self.species=='electron'.upper():
            self.mass_eV=const.physical_constants['electron mass energy equivalent in MeV'][0]*1e6
            self.charge=-1
        elif self.species=='proton'.upper():
            self.mass_eV=const.physical_constants['proton mass energy equivalent in MeV'][0]*1e6
            self.charge=1
        else:
            print('Unknown species {}'.format(self.species))
        self.classical_radius=1.0/4/np.pi/const.epsilon_0*const.e/self.mass_eV
        self.c_gamma=4.0*np.pi/3.0*self.classical_radius/np.power(self.mass_eV/1.0e6, 3.0)

    def get_gamma_beta(self, eng):
        gamma=eng/self.mass_eV
        beta=np.sqrt(1-1/gamma/gamma)
        return gamma,beta

    def get_energy(self, ref_eng, dpop):
        new_p=np.sqrt(ref_eng*ref_eng-self.mass_eV*self.mass_eV)*(1+dpop)
        new_e=np.sqrt(new_p*new_p+self.mass_eV*self.mass_eV)
        return new_e

    def get_momentum(self, eng):
        return np.sqrt(eng*eng-self.mass_eV*self.mass_eV)
