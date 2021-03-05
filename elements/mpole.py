from .. element import Element
from .. track import pass_components
import copy
import numpy as np

class Multipole(Element):
    SI_index_factor=0.175603595979828817023843904486
    SI_index_4th=[SI_index_factor+0.5, 2*SI_index_factor+1, -SI_index_factor,
                  -4*SI_index_factor-1, -SI_index_factor, 2*SI_index_factor+1,
                  SI_index_factor+0.5]
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
        self.properties['N_STEP'] = 100
        self.set_param(**param)

    def post_setting(self, name, value):
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
            self.properties['AERROR'][order] = fix_value + np.random.randn() * rms_value
        else:
            self.properties['BERROR'][order] = fix_value + np.random.randn() * rms_value
        pass

    def track(self, coor, math_lib):

        nstep = self.properties['N_STEP']
        #dis_sep = (int)(nstep / self.properties['N_OUT'])
        l = self.properties['L']
        an = self.properties['AN']
        bn = self.properties['BN']
        angle = self.properties['ANGLE']
        langle = self.properties['E1']
        rangle = self.properties['E2']

        SRout = 0
        if l == 0:
            if angle == 0:
                self.apply_align_error(coor, entrance=1.0)
                SRout = SRout + pass_components.pass_mpole(coor, 0, 0, an, bn, particle=None, ref_energy=0,
                                                      math_lib=math_lib, SR=0)
                self.apply_align_error(coor, entrance=-1.0)
                return SRout

            pass_components.pass_y_rot(coor, angle, math_lib=math_lib)

            return SRout

        cur = angle / l


        self.apply_align_error(coor, entrance=1.0)
        if len(bn) == 1:
            k1 = 0
        else:
            k1 = bn[1]

        pass_components.pass_dipole_fringe(coor, cur=cur, ang=langle)
        lstep=l/nstep
        for i in range(nstep):
            #ls += 1.0 * l / nstep
            for j in range(len(Multipole.SI_index_4th)):
                if j % 2 == 0:
                    pass_components.pass_drift(coor, Multipole.SI_index_4th[j] * lstep, math_lib=math_lib)
                else:
                    SRout = SRout + pass_components.pass_mpole(coor, Multipole.SI_index_4th[j] * lstep, cur, an, bn,
                                                          particle=None, ref_energy=0, math_lib=math_lib,
                                                          SR=0)


        pass_components.pass_dipole_fringe(coor, cur=cur, ang=rangle)
        self.apply_align_error(coor, entrance=-1.0)

        return SRout

'''
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
