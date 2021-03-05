from .. element import Element
from .mpole import Multipole
import copy

class Quadrupole(Multipole):
    propertyNames = copy.deepcopy(Multipole.propertyNames)

    # propertyNames.update(['L_ANGLE', 'R_ANGLE'])


    def __init__(self, name, **param):  # set type, defaults
        Multipole.__init__(self, name, DESIGN_ORDER=1)
        self.properties['K1'] = 0
        self.properties['BN'] = [0.0, 0.0]
        # self.properties['R_ANGLE']=0
        # self.properties['L_ANGLE']=0

        self.set_param(**param)


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