from .. element import Element
import copy

class Center(Element):
    propertyNames = copy.deepcopy(Element.propertyNames)
    propertyNames.update([
        'XC', 'PXC', 'YC', 'PYC', 'DEC', 'CTC'])

    def __init__(self, name, **param):  # set type, defaults
        Element.__init__(self, name, self.__class__.__name__)
        self.properties['XC'] = 0
        self.properties['PXC'] = 0
        self.properties['YC'] = 0
        self.properties['PYC'] = 0
        self.properties['DEC'] = 0
        self.properties['CTC'] = 0
        self.set_param(**param)