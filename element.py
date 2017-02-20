import copy
import numpy.random as random
import numpy as np


class Element(object):
    elementTypes = set(
        ['DRIFT', 'DIPOLE', 'QUADRUPOLE', 'SEXTRUPOLE', 'OCTUPOLE', 'MULTIPOLE', 'SOLENOID', 'MATRIX', 'TPSMAP',
         'BEAMBEAM', 'CAVITY', 'MARKER', 'BPM', 'SOLENOID', 'KICKER', 'CENTER', 'MALIGN', 'PROFILE','STRIPPER']
    )
    propertyNames = set(
        ['NAME', 'L', 'TYPE', 'RMS_DISPLACEMENT', 'RMS_ANGLE_DISPLACEMENT', 'TILT', 'N_STEP', 'APERTURE',
         'DISPLACEMENT', 'ANGLE_DISPLACEMENT', 'N_OUT']
    )

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
        self.properties['N_OUT'] = 1

        self.occurs = {}

    def getP(self, pname):
        if pname.upper() in self.properties:
            return self.properties[pname.upper()]
        else:
            print("The property name '{}' is not defined.".format(pname.upper()))

    def setP(self, **param):
        for k, v in param.items():
            if k.upper() in self.__class__.propertyNames:
                self.properties[k.upper()] = v
                self.postSetting(k.upper(), v)
            else:
                print('Uncognized property name {} in {} element'.format(k, self.__class__.__name__))
        self.check_consistency()

    def inBeamline(self, linename, pos):
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

    def getReverse(self, newname):
        return self.clone(newname)

    def check_consistency(self):
        if self.properties['N_OUT'] > self.properties['N_STEP']:
            self.properties['N_OUT'] = self.properties['N_STEP']

    def postSetting(self, name, value):
        pass




class Drift(Element):
    propertyNames = copy.deepcopy(Element.propertyNames)
    propertyNames.update(['EXACT_DRIFT'])
    def __init__(self, name, **param):
        Element.__init__(self, name, self.__class__.__name__)
        self.properties['EXACT_DRIFT'] = True
        self.setP(**param)

        # def build_CL_program(self, ctx):
        #    Drift.pass_cl_kernels= cl.Program(ctx, CL_track.kernel_pass_drift).build()

        # def track_GPU(self, coor_buf, ctx, queue, num_p, dim):
        #    Drift.pass_cl_kernels.pass_drift.set_scalar_arg_dtypes([None, np.float, np.uint32,np.uint32])
        #    Drift.pass_cl_kernels.pass_drift(coor_buf, self.properties['L'], num_p, dim)





class Mpole(Element):

    propertyNames = copy.deepcopy(Element.propertyNames)
    propertyNames.update([
        'BN', 'AN', 'BERROR', 'AERROR',
        'DESIGN_ORDER', 'MAX_ORDER', 'ANGLE', 'E1', 'E2',
        'K1', 'K2', 'K3', 'K4', 'K5', 'K6', 'K7', 'K8',
        'K1S', 'K2S', 'K3S', 'K4S', 'K5S', 'K6S', 'K7S', 'K8S'])

    def __init__(self, name, **param):  # set type, defaults
        Element.__init__(self, name, self.__class__.__name__)
        self.properties['BN'] = [0.0, ]
        self.properties['AN'] = [0.0, ]
        self.properties['BERROR'] = [0.0, ]
        self.properties['AERROR'] = [0.0, ]
        self.properties['DESIGN_ORDER'] = 0
        self.properties['ANGLE'] = 0
        self.properties['MAX_ORDER'] = 0
        self.properties['E1'] = 0
        self.properties['E2'] = 0
        self.setP(**param)

    def postSetting(self, name, value):
        if name[0] == 'K' and len(name) == 2:
            poleorder = int(name[1])
            for listname in ['BN', 'AN', 'BERROR', 'AERROR']:
                curlen = len(self.properties[listname])
                if curlen <= poleorder:
                    self.properties[listname].extend([0.0 for i in range(poleorder - curlen + 1)])
            self.properties['BN'][poleorder] = value
            if poleorder > self.properties['MAX_ORDER']:
                self.properties['MAX_ORDER'] = poleorder
        if name[0] == 'K' and len(name) == 3 and name[-1] == 'S':
            poleorder = int(name[1])
            for listname in ['BN', 'AN', 'BERROR', 'AERROR']:
                curlen = len(self.properties[listname])
                if curlen <= poleorder:
                    self.properties[listname].extend([0.0 for i in range(poleorder - curlen + 1)])
            self.properties['AN'][poleorder] = value
            if poleorder > self.properties['MAX_ORDER']:
                self.properties['MAX_ORDER'] = poleorder

        return

    def setMagnetError(self, order=0, rms_value=0, fix_value=0, skew=False):
        curlen = len(self.properties['BERROR'])
        if order >= curlen:
            self.properties['BERROR'].extend([0 for i in range(order - curlen + 1)])
            self.properties['AERROR'].extend([0 for i in range(order - curlen + 1)])
        if skew:
            self.properties['AERROR'][order] = fix_value + random.randn() * rms_value
        else:
            self.properties['BERROR'][order] = fix_value + random.randn() * rms_value
        pass

    '''def track(self, coor, particle, ref_energy, math_lib, ls, s_output=None, tps_list=None, centroid=None, SR=1):

        nstep = self.properties['N_STEP']
        dis_sep = (int)(nstep / self.properties['N_OUT'])
        l = self.properties['L']
        an = self.properties['AN']
        bn = self.properties['BN']
        angle = self.properties['ANGLE']
        langle = self.properties['E1']
        rangle = self.properties['E2']

        SRout = 0
        if l == 0:
            if angle == 0:
                self.apply_Error(coor, entrance=1.0)
                SRout = SRout + basictrack.pass_mpole(coor, 0, 0, an, bn, particle=particle, ref_energy=ref_energy,
                                                      math_lib=math_lib, SR=SR)
                self.apply_Error(coor, entrance=-1.0)
                if math_lib == np:
                    return SRout
                else:
                    t1, t2 = get_tpsa_linear_map(coor)
                    # t1,t2=coor.get_linear_map()
                    tps_list[-1] = (t1)
                    centroid[-1] = (t2)

                    return SRout
            basictrack.pass_y_rot(coor, angle, math_lib=math_lib)

            if math_lib == np:
                return SRout
            else:

                t1, t2 = get_tpsa_linear_map(coor)
                # t1,t2=coor.get_linear_map()
                tps_list[-1] = (t1)
                centroid[-1] = (t2)

                return SRout

        cur = angle / l
        if math_lib == np:
            pass
        else:
            if ls > s_output[-1]:
                s_output.append(ls)
                t1, t2 = get_tpsa_linear_map(coor)
                # t1,t2=coor.get_linear_map()
                tps_list.append(t1)
                centroid.append(t2)

        self.apply_Error(coor, entrance=1.0)
        if len(bn) == 1:
            k1 = 0
        else:
            k1 = bn[1]

        basictrack.pass_dipole_fringe(coor, cur=cur, ang=langle, math_lib=math_lib)

        for i in range(nstep):
            ls += 1.0 * l / nstep
            for j in range(len(Mpole.SI_index_4th)):
                if j % 2 == 0:
                    basictrack.pass_drift(coor, Mpole.SI_index_4th[j] * l / nstep, math_lib=math_lib)
                else:
                    SRout = SRout + basictrack.pass_mpole(coor, Mpole.SI_index_4th[j] * l / nstep, cur, an, bn,
                                                          particle=particle, ref_energy=ref_energy, math_lib=math_lib,
                                                          SR=SR)
            if (i + 1) % dis_sep == 0 or i == nstep - 1:

                if math_lib == np:
                    pass
                else:
                    s_output.append(ls)
                    self.apply_Error(coor, entrance=-1.0)
                    t1, t2 = get_tpsa_linear_map(coor)
                    # t1,t2=coor.get_linear_map()
                    self.apply_Error(coor, entrance=1.0)
                    tps_list.append(t1)
                    centroid.append(t2)

        basictrack.pass_dipole_fringe(coor, cur=cur, ang=rangle, math_lib=math_lib)
        self.apply_Error(coor, entrance=-1.0)

        return SRout

    @classmethod
    def build_CL_program(self, ctx):
        Mpole.pass_cl_kernels = cl.Program(ctx, cl_track.kernel_pass_drift + cl_track.kernel_pass_mpole).build()
        Mpole.pass_drift = Mpole.pass_cl_kernels.pass_drift
        Mpole.pass_drift.set_scalar_arg_dtypes([None, np.float, np.uint32, np.uint32])
        Mpole.pass_mpole = Mpole.pass_cl_kernels.pass_mpole
        Mpole.pass_mpole.set_scalar_arg_dtypes(
            [None, None, None, np.uint32, np.float, np.float, np.uint32, np.uint32, np.float, np.float])

    def track_cl(self, coor_buf, ctx, queue, num_p, dim, particle, ref_energy):

        nstep = self.properties['N_STEP']
        l = self.properties['L']
        an = self.properties['AN']
        bn = self.properties['BN']
        angle = self.properties['ANGLE']
        an_buf = cl.Buffer(ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(an))
        bn_buf = cl.Buffer(ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(bn))

        if l == 0:
            if angle == 0:
                Mpole.pass_mpole(queue, (num_p,), None, coor_buf, bn_buf, an_buf, np.int32(len(bn)), np.float(0.0),
                                 np.float(0.0), num_p, dim, np.float(ref_energy), np.float(particle.c_gamma))
                return

            return
        cur = angle / l
        for i in range(nstep):
            for j in range(len(Mpole.SI_index_4th)):
                if j % 2 == 0:

                    Mpole.pass_drift(queue, (num_p,), None, coor_buf, np.float(Mpole.SI_index_4th[j] * l / nstep),
                                     np.int32(num_p), np.int32(dim))
                else:
                    Mpole.pass_mpole(queue, (num_p,), None, coor_buf, bn_buf, an_buf, np.int32(len(bn)),
                                     np.float(Mpole.SI_index_4th[j] * l / nstep),
                                     np.float(cur), np.int32(num_p), np.int32(dim), np.float(ref_energy),
                                     np.float(particle.c_gamma))

'''


class Dipole(Mpole):
    propertyNames = copy.deepcopy(Mpole.propertyNames)
    propertyNames.update(['REF_ENERGY', 'EXACT_METHOD'])

    def __init__(self, name, **param):  # set type, defaults
        Mpole.__init__(self, name, DESIGN_ORDER=0)
        self.properties['E1'] = 0
        self.properties['E2'] = 0
        self.setP(**param)

    def getReverse(self, newname):
        temp = self.clone(newname)
        if 'E1' in self.properties:
            E1copy = self.getP('E1')
            temp.setP(E1=self.getP('E2'))
            temp.setP(E2=E1copy)
        if 'E2' in self.properties:
            E1copy = self.getP('E1')
            temp.setP(E1=self.getP('E2'))
            temp.setP(E2=E1copy)
        return temp


class Quadrupole(Mpole):
    propertyNames = copy.deepcopy(Mpole.propertyNames)

    # propertyNames.update(['L_ANGLE', 'R_ANGLE'])


    def __init__(self, name, **param):  # set type, defaults
        Mpole.__init__(self, name, DESIGN_ORDER=1)
        self.properties['K1'] = 0
        self.properties['BN'] = [0.0, 0.0]
        # self.properties['R_ANGLE']=0
        # self.properties['L_ANGLE']=0

        self.setP(**param)

    def check_consistency(self):
        Element.check_consistency(self)
        # if self.properties['R_ANGLE']+self.properties['L_ANGLE']==self.properties['ANGLE']:
        #    pass
        # else:
        #    print("Something Wring with LR angle and angle")
        #    exit()

'''
    def track(self, coor, particle, ref_energy, math_lib, ls, s_output=None, tps_list=None, centroid=None, SR=1):

        langle = self.properties['E1']
        rangle = self.properties['E2']

        nstep = self.properties['N_STEP']
        dis_sep = (int)(nstep / self.properties['N_OUT'])
        l = self.properties['L']
        an = self.properties['AN']
        bn = self.properties['BN']
        SRout = 0

        if math_lib == np:
            pass
        else:
            if ls > s_output[-1]:
                s_output.append(ls)
                t1, t2 = get_tpsa_linear_map(coor)
                # t1,t2=coor.get_linear_map()
                tps_list.append(t1)
                centroid.append(t2)

        basictrack.pass_y_rot(coor, langle, math_lib=math_lib)

        self.apply_Error(coor, entrance=1.0)
        basictrack.pass_quad_fringe(coor, bn[1], math_lib=math_lib)
        if l == 0:
            SRout = SRout + basictrack.pass_mpole(coor, 0, 0, an, bn, particle=particle, ref_energy=ref_energy,
                                                  math_lib=math_lib, SR=SR)
            if math_lib == np:
                pass
            else:
                self.apply_Error(coor, entrance=-1.0)
                t1, t2 = get_tpsa_linear_map(coor)
                # t1,t2=coor.get_linear_map()
                self.apply_Error(coor, entrance=1.0)
                tps_list[-1] = t1
                centroid[-1] = t2
        else:
            for i in range(nstep):
                ls += 1.0 * l / nstep
                for j in range(len(Mpole.SI_index_4th)):
                    if j % 2 == 0:
                        basictrack.pass_drift(coor, Mpole.SI_index_4th[j] * l / nstep, math_lib=math_lib)
                    else:
                        SRout = SRout + basictrack.pass_mpole(coor, Mpole.SI_index_4th[j] * l / nstep, 0, an, bn,
                                                              particle=particle, ref_energy=ref_energy,
                                                              math_lib=math_lib, SR=SR)

                if (i + 1) % dis_sep == 0 and i < nstep - 1:
                    if math_lib == np:
                        pass
                    else:
                        s_output.append(ls)
                        self.apply_Error(coor, entrance=-1.0)
                        # t1,t2=coor.get_linear_map()
                        t1, t2 = get_tpsa_linear_map(coor)
                        self.apply_Error(coor, entrance=1.0)
                        tps_list.append(t1)
                        centroid.append(t2)
        basictrack.pass_quad_fringe(coor, -bn[1], math_lib=math_lib)
        self.apply_Error(coor, entrance=-1.0)

        basictrack.pass_y_rot(coor, rangle, math_lib=math_lib)
        if math_lib == np:
            pass
        else:
            s_output.append(ls)
            t1, t2 = get_tpsa_linear_map(coor)

            # t1,t2=coor.get_linear_map()
            tps_list.append(t1)
            centroid.append(t2)

        return SRout

    @classmethod
    def build_CL_program(self, ctx):
        Quadrupole.pass_cl_kernels = cl.Program(ctx,
                                                cl_track.kernel_pass_drift + cl_track.kernel_pass_mpole + cl_track.kernel_pass_y_rot).build()
        Quadrupole.pass_drift = Quadrupole.pass_cl_kernels.pass_drift
        Quadrupole.pass_drift.set_scalar_arg_dtypes([None, np.float, np.uint32, np.uint32])
        Quadrupole.pass_mpole = Quadrupole.pass_cl_kernels.pass_mpole
        Quadrupole.pass_mpole.set_scalar_arg_dtypes(
            [None, None, None, np.uint32, np.float, np.float, np.uint32, np.uint32, np.float, np.float])
        Quadrupole.pass_y_rot = Quadrupole.pass_cl_kernels.pass_y_rot
        Quadrupole.pass_y_rot.set_scalar_arg_dtypes([None, np.float, np.uint32, np.uint32])

    def track_cl(self, coor_buf, ctx, queue, num_p, dim, particle, ref_energy):
        nstep = self.properties['N_STEP']
        l = self.properties['L']
        an = self.properties['AN']
        bn = self.properties['BN']
        angle = self.properties['ANGLE']
        an_buf = cl.Buffer(ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(an))
        bn_buf = cl.Buffer(ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(bn))

        if angle != 0:
            Quadrupole.pass_y_rot(queue, (num_p,), None, coor_buf, np.float(angle / 2.0), num_p, dim)
        if l == 0:
            Quadrupole.pass_mpole(queue, (num_p,), None, coor_buf, bn_buf, an_buf, np.uint32(len(bn)), np.float(0.0),
                                  np.float(0.0), num_p, dim, np.float(ref_energy), np.float(particle.c_gamma))
        else:
            for i in range(nstep):
                for j in range(len(Mpole.SI_index_4th)):
                    if j % 2 == 0:

                        Quadrupole.pass_drift(queue, (num_p,), None, coor_buf,
                                              np.float(Mpole.SI_index_4th[j] * l / nstep), num_p, dim)
                    else:
                        Quadrupole.pass_mpole(queue, (num_p,), None, coor_buf, bn_buf, an_buf, np.uint32(len(bn)),
                                              np.float(Mpole.SI_index_4th[j] * l / nstep),
                                              np.float(0.0), num_p, dim, np.float(ref_energy),
                                              np.float(particle.c_gamma))
        if angle != 0:
            Quadrupole.pass_y_rot(queue, (num_p,), None, coor_buf, np.float(angle / 2.0), num_p, dim)
'''

class Sextrupole(Mpole):
    propertyNames = copy.deepcopy(Mpole.propertyNames)

    def __init__(self, name, **param):  # set type, defaults
        Mpole.__init__(self, name, DESIGN_ORDER=2)
        self.properties['K2'] = 0
        self.properties['BN'] = [0.0, 0.0, 0.0]

        self.setP(**param)


class Octupole(Mpole):
    propertyNames = copy.deepcopy(Mpole.propertyNames)

    def __init__(self, name, **param):  # set type, defaults
        Mpole.__init__(self, name, DESIGN_ORDER=3)
        self.properties['K3'] = 0
        self.properties['BN'] = [0.0, 0.0, 0.0, 0.0]

        self.setP(**param)


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
        self.setP(**param)
'''
    def track(self, coor, particle, ref_energy, math_lib, ls, s_output, tps_list, SR=1):
        pass
'''

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
        self.setP(**param)


class MALIGN(Element):
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


class Bpm(Element):
    propertyNames = copy.deepcopy(Element.propertyNames)
    propertyNames.update([
        'WEIGHT', ])

    def __init__(self, name, **param):  # set type, defaults
        Element.__init__(self, name, self.__class__.__name__)
        self.properties['WEIGHT'] = 1
        self.setP(**param)


class Marker(Element):
    propertyNames = copy.deepcopy(Element.propertyNames)

    def __init__(self, name, **param):  # set type, defaults
        Element.__init__(self, name, self.__class__.__name__)
        self.setP(**param)


class Kicker(Element):
    propertyNames = copy.deepcopy(Element.propertyNames)
    propertyNames.update([
        'DPX', 'DPY'])

    def __init__(self, name, **param):  # set type, defaults
        Element.__init__(self, name, self.__class__.__name__)
        self.properties['DPX'] = 0.0
        self.properties['DPY'] = 0.0
        self.setP(**param)


class Matrix(Element):
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



