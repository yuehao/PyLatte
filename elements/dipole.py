from .mpole import Multipole
import copy

class Dipole(Multipole):
    propertyNames = copy.deepcopy(Multipole.propertyNames)
    propertyNames.update(['REF_ENERGY', 'EXACT_METHOD'])

    def __init__(self, name, **param):  # set type, defaults
        Multipole.__init__(self, name, DESIGN_ORDER=0)
        self.properties['E1'] = 0
        self.properties['E2'] = 0
        self.set_param(**param)

    def getReverse(self, newname):
        temp = self.clone(newname)
        if 'E1' in self.properties:
            E1copy = self.get_param('E1')
            temp.set_param(E1=self.get_param('E2'))
            temp.set_param(E2=E1copy)
        if 'E2' in self.properties:
            E1copy = self.get_param('E1')
            temp.set_param(E1=self.get_param('E2'))
            temp.set_param(E2=E1copy)
        return temp