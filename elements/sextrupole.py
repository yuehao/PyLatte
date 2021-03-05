from .mpole import Multipole
import copy

class Sextrupole(Multipole):
    propertyNames = copy.deepcopy(Multipole.propertyNames)

    def __init__(self, name, **param):  # set type, defaults
        Multipole.__init__(self, name, DESIGN_ORDER=2)
        self.properties['K2'] = 0
        self.properties['BN'] = [0.0, 0.0, 0.0]

        self.set_param(**param)