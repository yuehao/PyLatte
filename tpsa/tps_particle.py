import PyTPSA as ptps
import numpy as np

from ..particle import Particle


class tps_particle(object):
    initilized = False
    dimension = 0
    max_order = 0
    @classmethod
    def initilize(cls, dim, order):
        ptps.tpsa.initilize(dim, order)
        tps_particle.initilized = True
        tps_particle.dimension = dim
        tps_particle.max_order = order
    def __init__(self, ini_value = np.zeros(6),
                 default_dim=6, default_order=3):
        if tps_particle.initilized is not True:
            tps_particle.initilize(default_dim, default_order)

        self.dimension=tps_particle.dimension

        self.ini_value = ini_value
        self.particle = Particle(species='electron')

        self.coordinates = []
        for i in range(self.dimension):
            self.coordinates.append(ptps.tpsa(self.ini_value[i],i+1))

    def get_value(self):
        self.value=np.array([t.cst() for t in self.coordinates])
        return self.value

    def get_linear_map(self, scalar=False):
        if scalar:
            linear_map=[self.coordinates[i].element(j+1) for i in range(self.dimension) for j in range(self.dimension)]
        else:
            linear_map=[self.coordinates[i].derivative(j+1) for i in range(self.dimension) for j in range(self.dimension)]

        #consts=np.array([self.coordinates[i].cst() for i in range(self.dimension)])
        return np.array(linear_map).reshape(self.dimension,self.dimension)
