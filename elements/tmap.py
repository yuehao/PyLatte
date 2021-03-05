from .. element import Element
import copy
import numpy as np

class Tmap(Element):
    propertyNames = copy.deepcopy(Element.propertyNames)
    propertyNames.update([
        'D1', 'D2', 'D3', 'D4', 'D5', 'D6',
        'R11', 'R12', 'R13', 'R14', 'R15', 'R16',
        'R21', 'R22', 'R23', 'R24', 'R25', 'R26',
        'R31', 'R32', 'R33', 'R34', 'R35', 'R36',
        'R41', 'R42', 'R43', 'R44', 'R45', 'R46',
        'R51', 'R52', 'R53', 'R54', 'R55', 'R56',
        'R61', 'R62', 'R63', 'R64', 'R65', 'R66'])

    def __init__(self, name, dim=6, **param):
        Element.__init__(self, name, self.__class__.__name__)
        self.dimension = dim
        self.vector = np.zeros(dim)
        self.matrix = np.eye(dim)
        for i in range(6):
            cname = 'D{}'.format(i)
            self.properties[cname] = 0.0
            for j in range(6):
                matname = 'R{}{}'.format(i, j)
                if i == j:
                    self.properties[matname] = 1.0
                else:
                    self.properties[matname] = 0.0

    def postSetting(self, name, value):
        if name[0] == 'D' and len(name) == 2 and int(name[1]) <= self.dimension:
            self.vector[int(name[1]) - 1] = value
        if name[0] == 'R' and len(name) == 3 and int(name[1]) <= self.dimension and int(name[2]) <= self.dimension:
            self.vector[int(name[1]) - 1, int(name[2]) - 1] = value

