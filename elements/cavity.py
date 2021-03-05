from .. element import Element
import copy


class Cavity(Element):
    propertyNames = copy.deepcopy(Element.propertyNames)
    propertyNames.update([
        'VOLT', 'FREQ', 'PHASE', 'PHASE_REFERENCE', 'PHASE_ADJUSTMENT', 'HOM_FREQ', 'HOM_Q', 'HOM_ROQ', 'LINAC_MODEL',
        'UPDATE_ENERGY'])
    timeDependent = True

    def __init__(self, name, **param):  # set type, defaults
        Element.__init__(self, name, self.__class__.__name__)
        self.properties['VOLT'] = 0
        self.properties['FREQ'] = 1e9
        self.properties['PHASE'] = 0
        self.properties['UPDATE_ENERGY'] = 1
        self.properties['PHASE_REFERENCE'] = 0
        self.properties['PHASE_ADJUSTMENT'] = 0
        self.properties['HOM_FREQ'] = []
        self.properties['HOM_Q'] = []
        self.properties['HOM_ROQ'] = []
        self.properties['LINAC_MODEL'] = 'SRS'
        self.set_param(**param)
'''
    def track(self, coor, particle, ref_energy, math_lib, ls, s_output, tps_list, SR=1):
        pass
'''
