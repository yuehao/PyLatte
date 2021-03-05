from .. element import Element
import copy


class Marker(Element):
    propertyNames = copy.deepcopy(Element.propertyNames)

    def __init__(self, name, **param):  # set type, defaults
        Element.__init__(self, name, self.__class__.__name__)
        self.setP(**param)