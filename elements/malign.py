from .. element import Element
import copy


class Malign(Element):
    propertyNames = copy.deepcopy(Element.propertyNames)
    propertyNames.update([
        'DX', 'DPX', 'DY', 'DPY', 'DEOE', 'DCT'])

    def __init__(self, name, **param):  # set type, defaults
        Element.__init__(self, name, self.__class__.__name__)
        self.properties['DX'] = 0
        self.properties['DPX'] = 0
        self.properties['DY'] = 0
        self.properties['DPY'] = 0
        self.properties['DEOE'] = 0
        self.properties['DCT'] = 0
        self.setP(**param)


