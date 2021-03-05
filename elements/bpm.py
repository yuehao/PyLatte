from .. element import Element
import copy


class Bpm(Element):
    propertyNames = copy.deepcopy(Element.propertyNames)
    propertyNames.update([
        'WEIGHT', ])

    def __init__(self, name, **param):  # set type, defaults
        Element.__init__(self, name, self.__class__.__name__)
        self.properties['WEIGHT'] = 1
        self.set_param(**param)