from .mpole import Multipole
import copy



class Octupole(Multipole):
    propertyNames = copy.deepcopy(Multipole.propertyNames)

    def __init__(self, name, **param):  # set type, defaults
        Multipole.__init__(self, name, DESIGN_ORDER=3)
        self.properties['K3'] = 0
        self.properties['BN'] = [0.0, 0.0, 0.0, 0.0]

        self.set_param(**param)