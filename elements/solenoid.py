from .. element import Element
import copy

class Solenoid(Element):
    propertyNames = copy.deepcopy(Element.propertyNames)
    propertyNames.update([
        'B', 'KS'])

    def __init__(self, name, **param):  # set type, defaults
        Element.__init__(self, name, self.__class__.__name__)
        self.properties['B'] = 0.0
        self.properties['KS'] = 0.0
        self.set_param(**param)