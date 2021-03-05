import numpy as np
from .particle import Particle
from .optics import Optics

def random_Gaussian(dim, cutoff=4, seed=None):
    np.random.seed(seed)
    def get_One(cutoff):
        while True:
            res=np.random.randn()
            if abs(res)<cutoff:
                yield res

    iter=get_One(cutoff)
    dist=np.array([iter.next() for i in range(dim)])

    average=np.average(dist)
    dist-=average
    size=np.std(dist)
    if size>0:
        dist/=size

    return dist


_x_, _px_, _y_, _py_,_z_,_de_=range(6)

class Bunch(object):
    def __init__(self, species='electron', **param):

        self.num_particle=1
        self.dimension=6
        self.particle=Particle(species)
        self.bunch_spacing=1e-6  #default bunch spacing 1 us
        self.ini_energy=1.0e9

        self.ini_momentum=np.sqrt(self.ini_energy*self.ini_energy-self.particle.mass_eV*self.particle.mass_eV)
        self.ini_emittance_x=1e-6
        self.ini_emittance_y=1e-6
        self.ini_bunch_length=1e-3
        self.ini_energy_spread=1e-5
        self.ini_optics=Optics()
        self.coordinates=np.zeros((self.dimension, self.num_particle))
        self.setP(**param)


    def setP(self, **param):
        for k,v in param.iteritems():
            if k.lower() in self.__dict__:
                self.__dict__[k.lower()] = v
                if k.lower()=='num_particle' or k.lower()=='dimension':
                    self.coordinates=np.zeros((self.dimension, self.num_particle))
            else:
                print('Unrecognized bunch parameter {} with its value {}'.format(k.lower(), v))

    def set_ini_distribution(self, cutoff=4, reference_coor=None):

        self.coordinates[_x_,:]=random_Gaussian(self.num_particle)*np.sqrt(self.ini_emittance_x*self.ini_optics.beta_x)
        self.coordinates[_px_, :]=random_Gaussian(self.num_particle)*np.sqrt(self.ini_emittance_x*self.ini_optics.beta_x)/self.ini_optics.beta_x
        self.coordinates[_px_, :]+=self.coordinates[_x_, :]*self.ini_optics.alpha_x/self.ini_optics.beta_x
        self.coordinates[_y_, :]=random_Gaussian(self.num_particle)*np.sqrt(self.ini_emittance_y*self.ini_optics.beta_y)
        self.coordinates[_py_, :]=random_Gaussian(self.num_particle)*np.sqrt(self.ini_emittance_y*self.ini_optics.beta_y)/self.ini_optics.beta_y
        self.coordinates[_py_, :]+=self.coordinates[_y_, :]*self.ini_optics.alpha_y/self.ini_optics.beta_y

        self.coordinates[_z_, :]=random_Gaussian(self.num_particle)*self.ini_bunch_length*self.particle.get_gamma_beta(self.ini_energy)[1]
        self.coordinates[_de_, :]=random_Gaussian(self.num_particle)*self.ini_energy_spread

        self.coordinates[_x_, :]+=self.coordinates[_de_, :]*self.ini_optics.eta_x
        self.coordinates[_px_, :]+=self.coordinates[_de_, :]*self.ini_optics.eta_xp
        self.coordinates[_y_, :]+=self.coordinates[_de_, :]*self.ini_optics.eta_y
        self.coordinates[_py_, :]+=self.coordinates[_de_, :]*self.ini_optics.eta_yp

        if reference_coor is not None and len(reference_coor)==self.dimension:

            self.coordinates+=np.outer(reference_coor, np.ones(self.num_particle))


    def statistics(self, ref_energy):
        centroid=np.mean(self.coordinates, axis=1)
        moments=np.zeros((self.dimension,self.dimension))

        for i in range(self.dimension):
            for j in range(i, self.dimension):
                moments[i,j]=np.correlate(self.coordinates[i,:]-centroid[i],self.coordinates[j,:]-centroid[j])/self.num_particle
        for i in range(1,self.dimension):
                for j in range(i):
                    moments[i,j]=moments[j,i]


        average_p=np.sqrt(ref_energy*ref_energy-self.particle.mass_eV*self.particle.mass_eV)*(1+centroid[_de_])
        average_e=np.sqrt(average_p*average_p+self.particle.mass_eV*self.particle.mass_eV)
        ave_gamma_beta=average_p/self.particle.mass_eV

        emit_x=np.sqrt(moments[_x_,_x_]*moments[_px_,_px_]-moments[_x_,_px_]*moments[_x_,_px_])
        norm_emit_x=emit_x*ave_gamma_beta
        emit_y=np.sqrt(moments[_y_,_y_]*moments[_py_,_py_]-moments[_y_,_py_]*moments[_y_,_py_])
        norm_emit_y=emit_y*ave_gamma_beta

        return centroid, moments, average_e, emit_x, norm_emit_x, emit_y, norm_emit_y

