from .. element import Element
import copy


class Kicker(Element):
    propertyNames = copy.deepcopy(Element.propertyNames)
    propertyNames.update([
        'DPX', 'DPY'])

    def __init__(self, name, **param):  # set type, defaults
        Element.__init__(self, name, self.__class__.__name__)
        self.properties['DPX'] = 0.0
        self.properties['DPY'] = 0.0
        self.set_param(**param)