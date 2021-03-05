import copy
import numpy.random as random
import numpy as np
from .track import pass_components


class Element(object):
    elementTypes = ['DRIFT', 'DIPOLE', 'QUADRUPOLE', 'SEXTRUPOLE', 'OCTUPOLE', 'MULTIPOLE', 'SOLENOID', 'TMAP',
         'CAVITY', 'MARKER', 'BPM', 'KICKER', 'CENTER', 'MALIGN']

    propertyNames = set(
        ['NAME', 'L', 'TYPE', 'RMS_DISPLACEMENT', 'RMS_ANGLE_DISPLACEMENT', 'TILT', 'N_STEP', 'APERTURE',
         'DISPLACEMENT', 'ANGLE_DISPLACEMENT', 'N_SPLIT']
    )

    @staticmethod
    def get_property_name(eletype):
        import importlib
        temp=importlib.import_module('.elements.'+eletype.lower(), package='PyLatte').__getattribute__(eletype.title())
        return temp.propertyNames


    def __init__(self, name, eletype):
        if eletype.upper() not in Element.elementTypes:
            print("The element has unrecognized type {}".format(eletype))
            return
        self.properties = {}
        self.properties['TYPE'] = eletype.upper()
        self.properties['L'] = 0
        self.properties['N_STEP'] = 1
        self.properties['APERTURE'] = 1
        self.properties['NAME'] = name.upper()
        self.properties['RMS_DISPLACEMENT'] = [0, 0, 0]
        self.properties['RMS_ANGLE_DISPLACEMENT'] = [0, 0, 0]
        self.properties['DISPLACEMENT'] = [0, 0, 0]
        self.properties['ANGLE_DISPLACEMENT'] = [0, 0, 0]  # the angles are theta(0,pi), phi(-pi,pi) around 0 and tilt error [-pi,pi] around tilt angle
        self.properties['TILT'] = 0
        self.properties['N_OUT'] =1


        self.occurs = {}

    def get_param(self, pname):
        '''

        :param pname: Property names
        :return: Property value
        '''
        if pname.upper() in self.properties:
            return self.properties[pname.upper()]
        else:
            print("The property name '{}' is not defined.".format(pname.upper()))

    def set_param(self, **param):
        '''

        :param param:
        :return:
        '''
        for k, v in param.items():
            if k.upper() == 'TYPE':
                #Does not allow Type to be redefined
                pass
            elif k.upper() in self.__class__.propertyNames:
                self.properties[k.upper()] = v
                self.post_setting(k.upper(), v)
            else:
                print('Uncognized property name {} in {} element'.format(k, self.__class__.__name__))
        self.check_consistency()
    def track(self,coor, math_lib):
        pass
    def in_beamline(self, linename, pos):
        if linename in self.occurs:
            self.occurs[linename].append(pos)
        else:
            self.occurs[linename] = []
            self.occurs[linename].append(pos)

    def freeze_properties(self):
        temp = copy.deepcopy(self.properties)
        for k, v in temp.items():
            if isinstance(v, list):
                temp[k] = tuple(v)
        return temp

    def compare(self, other_ele):
        if len(self.properties) != len(other_ele.properties):
            return 0

        all_p = (set(self.freeze_properties().items()) and set(other_ele.freeze_properties().items()))
        if len(self.properties) == len(all_p):
            return 1
        return 0

    def clone(self, newname, include_error=False):
        temp = copy.deepcopy(self)
        temp.properties['NAME'] = newname
        temp.properties['GROUP'] = ''
        temp.occurs = {}
        # don't clone error
        if include_error:
            temp.properties['DISPLACEMENT'] = [0, 0, 0]
            temp.properties['ANGLE_DISPLACEMENT'] = [0, 0, 0]  # the angles are theta(0,pi), phi(-pi,pi) around 0 and tilt error [-pi,pi] around tilt angle
        return temp

    def get_reverse(self, newname):
        return self.clone(newname)

    def check_consistency(self):
        if self.properties['N_OUT'] > self.properties['N_STEP']:
            self.properties['N_OUT'] = self.properties['N_STEP']


    def post_setting(self, name, value):
        pass

    def apply_align_error(self, coor, entrance=1.0):
        if entrance == 1.0:
            pass_components.apply_transverse_offset(coor, self.properties['DISPLACEMENT'], apply_factor=entrance)
            if self.properties['ANGLE_DISPLACEMENT'][-1] != 0:
                pass_components.apply_rotation(coor, -self.properties['ANGLE_DISPLACEMENT'][-1], apply_factor=entrance)
        elif entrance == -1.0:
            if self.properties['ANGLE_DISPLACEMENT'][-1] != 0:
                pass_components.apply_rotation(coor, -self.properties['ANGLE_DISPLACEMENT'][-1], apply_factor=entrance)
            pass_components.apply_transverse_offset(coor, self.properties['DISPLACEMENT'], apply_factor=entrance)

















