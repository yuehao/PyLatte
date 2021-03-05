import numpy as np

_x_, _px_, _y_, _py_,_z_,_de_=range(6)
def sqr_ps2(coor):
    '''
    calculate the square of the ps
    :param coor: coordinate of the whole bunch coor[i,0] is the x of the ith particle
    :param exact_Hamiltonian: 1, Hamiltonian for [x,px,y,py,z=-\beta ct, delta=dp/p]

    :return:
    '''
    return (1.0+coor[_de_])*(1.0+coor[_de_])-coor[_px_]*coor[_px_]-coor[_py_]*coor[_py_]

def pass_drift(coor, l, Hamiltonian=1, math_lib=np):
    temp=l/math_lib.sqrt(sqr_ps2(coor))
    coor[_x_]+=coor[_px_]*temp
    coor[_y_]+=coor[_py_]*temp
    coor[_z_]+=((1+coor[_de_])*temp-l)

def pass_matrix(coor, mat):
    return np.dot(mat, coor)

def classical_radiation(coor, l, cur, bx, by, particle, ref_energy, Hamiltonian=1, math_lib=np):

    ps=math_lib.sqrt(sqr_ps2(coor))
    e_vec=np.array([coor[_px_],coor[_py_], (1.0+coor[_x_]*cur)*ps])

    abs_e_vec=math_lib.sqrt(coor[_px_]*coor[_px_]+coor[_py_]*coor[_py_]+(1.0+coor[_x_]*cur)*ps*(1.0+coor[_x_]*cur)*ps)
    e_vec=e_vec/abs_e_vec
    #b_vec=np.array([bx,by,0.0])

    b2=e_vec[2]*by*e_vec[2]*by+e_vec[2]*bx*e_vec[2]*bx+(e_vec[0]*by-e_vec[1]*bx)*(e_vec[0]*by-e_vec[1]*bx)


    energy_mev=ref_energy/1.0e6*(1+coor[_de_])
    deltaE=particle.c_gamma/2.0/np.pi*energy_mev*energy_mev*energy_mev*b2*(1+coor[_x_]*cur)*l/ps
    #print(deltaE.shape)
    coor[_de_]-=deltaE
    new_ps=math_lib.sqrt(sqr_ps2(coor))
    er=(new_ps/ps)

    coor[_px_]*=er
    coor[_py_]*=er
    return ref_energy*deltaE


def apply_transverse_offset(coor, error, apply_factor=1):
    if error[0] !=0:
        coor[_x_]-=(apply_factor*error[0])
    if error[1] !=0:
        coor[_y_]-=(apply_factor*error[1])



def apply_rotation(coor, ang, apply_factor=1):
    matrot=np.eye(4)
    angle=ang*apply_factor
    matrot[_x_, _x_]=np.cos(angle)
    matrot[_y_, _y_]=np.cos(angle)
    matrot[_x_, _y_]=-np.sin(angle)
    matrot[_y_, _x_]=np.sin(angle)

    matrot[_px_, _px_]=np.cos(angle)
    matrot[_py_, _py_]=np.cos(angle)
    matrot[_px_, _py_]=-np.sin(angle)
    matrot[_py_, _px_]=np.sin(angle)

    coor[0:4]=np.dot(matrot, coor[0:4])


def pass_y_rot(coor, ang, Hamiltonian=1, math_lib=np):

    from math import sin,cos,tan
    tan_ang=tan(ang)
    cos_ang=cos(ang)
    ps=math_lib.sqrt(sqr_ps2(coor))
    temp1=1.0/(ps-coor[_px_]*tan_ang)
    temp2=coor[_x_]*tan_ang
    coor[_px_]=coor[_px_]*cos_ang+sin(ang)*ps
    coor[_y_]+=coor[_py_]*temp2*temp1
    coor[_z_]+=(1.0+coor[_de_])*temp2*temp1
    coor[_x_]=coor[_x_]*ps*temp1/cos_ang

def pass_quad_fringe(coor, k1, cur=0, Hamiltonian=1, math_lib=np):
    ps=math_lib.sqrt(sqr_ps2(coor))
    ang=coor[_px_]/ps
    focusing_str=(k1*coor[_x_]+cur)*math_lib.tan(ang)
    coor[_px_]+=focusing_str*coor[_x_]
    coor[_py_]-=focusing_str*coor[_y_]


def pass_dipole_fringe(coor, cur, ang):
    from math import tan
    focusing_str=(cur)*tan(ang)
    coor[_px_]+=focusing_str*coor[_x_]
    coor[_py_]-=focusing_str*coor[_y_]

def pass_mpole(coor, l, cur, an, bn, particle=None, ang=0, ref_energy=0, Hamiltonian=1, math_lib=np, SR=0):
        '''
        :param coor: Coordinates
        :param l: length
        :param cur: curvature of the bending
        :param an: skew components, an = (d^nBx/dx^n)/Brho
        :param bn: normal components, bn = (d^nBy/dx^n)/Brho
        :param particle:
        :param ang:
        :param ref_energy:
        :param Hamiltonian:
        :param math_lib: numpy or PyTPSA
        :param SR:
        :return:
        '''
        import math
        if cur!=0:
            ps=math_lib.sqrt(sqr_ps2(coor))
            coor[_px_]+=(cur*(ps-1.0)-cur*cur*coor[_x_])*l
            coor[_z_]+=cur*l*coor[_x_]

        kreal=bn[0]
        kimag=an[0]
        accureal=1
        accuimag=0
        for i in range(1, len(bn)):
            accureal,accuimag=accureal*coor[_x_]-accuimag*coor[_y_], accuimag*coor[_x_]+accureal*coor[_y_]

            kreal+=(accureal*bn[i]-accuimag*an[i])/math.factorial(i)
            kimag+=(accureal*an[i]+accuimag*bn[i])/math.factorial(i)
        sr_energy=0
        if SR==1:
            sr_energy=classical_radiation(coor, l, cur, kimag, kreal+cur, particle=particle, ref_energy=ref_energy, Hamiltonian=Hamiltonian, math_lib=math_lib)
        coor[_px_]-=kreal*l
        coor[_py_]+=kimag*l

        return sr_energy

def pass_cavity_simple(coor, voltage, freq, phase):
    pass


def sqr_ps2_px2(coor):
    '''
    calculate the square of the ps
    :param coor: coordinate of the whole bunch coor[i,0] is the x of the ith particle
    :param exact_Hamiltonian: 1, Hamiltonian for [x,px,y,py,z=-\beta ct, delta=dp/p]

    :return:
    '''
    return (1.0+coor[_de_])*(1.0+coor[_de_])-coor[_py_]*coor[_py_]
def pass_dipole_geo(coor, cur, l, math_lib):
    psx2=sqr_ps2_px2(coor)
    from math import sin,cos
    sintheta=sin(cur*l)
    costheta=cos(cur*l)

    ps_old=math_lib.sqrt(psx2-coor[_px_]*coor[_px_])
    coor_2=coor[_px_]*costheta+(ps_old-1-coor[_x_]*cur)*sintheta
    coor_2_d = -cur * coor[_px_] * sintheta + cur * (ps_old - 1 - coor[_x_] * cur) * costheta
    ps = math_lib.sqrt(psx2 - coor_2 * coor_2)










