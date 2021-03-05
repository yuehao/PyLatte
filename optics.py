import numpy as np


x_, px_, y_, py_,z_,de_=range(6)

def total_derivative(y, x_ind, tps6d, exclude_list=[]):
    '''

    :param y: tps variable to be differentiated
    :param x_ind: the x index starting with zero
    :param tps6d: the full map
    :return:
    '''
    dim=len(tps6d)
    res=y.derivative(x_ind+1)
    for i in range(dim):
        if i!=x_ind and (i not in exclude_list):
            res+=y.derivative(i+1)*tps6d[i].derivative(x_ind+1)
    return res


class Optics(object):
    def __init__(self, math_library=np, **param):
        self.beta_x=1
        self.alpha_x=0
        self.gamma_x=1
        self.tune_x=0
        self.matched_tune_x=0

        self.beta_y=1
        self.alpha_y=0
        self.gamma_y=1
        self.tune_y=0
        self.matched_tune_y=0


        self.eta_x=0
        self.eta_xp=0
        self.eta_y=0
        self.eta_yp=0

        self.H_func_x=0
        self.H_func_y=0

        self.pass_diff=0
        self.c_alpha_c=0
        self.c_alpha_c_2=0
        self.chromaticity_x=0
        self.chromaticity_x_2=0
        self.chromaticity_y=0
        self.chromaticity_y_2=0

        self.setP(**param)
        self.norm_matrix(math_library=math_library)

    def setP(self, **param):
        for k,v in param.items():
            if k.lower() in self.__dict__:
                self.__dict__[k.lower()] = v
            else:
                print('Unrecognized optics name {} with its value {}'.format(k.lower(), v))
        self.gamma_x=(1+self.alpha_x*self.alpha_x)/self.beta_x
        self.gamma_y=(1+self.alpha_y*self.alpha_y)/self.beta_y
        self.H_func_x=self.gamma_x*self.eta_x*self.eta_x+2*self.alpha_x*self.eta_x*self.eta_xp+self.beta_x*self.eta_xp*self.eta_xp
        self.H_func_y=self.gamma_y*self.eta_y*self.eta_y+2*self.alpha_y*self.eta_y*self.eta_yp+self.beta_y*self.eta_yp*self.eta_yp


    def norm_matrix(self, math_library=np):
        self.norm_mat=np.eye(6)
        if math_library is not np:
            dfv=tpsa.CTPS(0.0)
            one=tpsa.CTPS(1.0)
            self.norm_mat=np.full((6,6), fill_value=dfv, dtype=object)
            self.norm_mat[z_,z_]=one
            self.norm_mat[de_,de_]=one

            #print(math_library.sqrt(self.beta_x))
        self.norm_mat[x_,x_]=math_library.sqrt(self.beta_x)
        self.norm_mat[px_,x_]=-self.alpha_x/math_library.sqrt(self.beta_x)
        self.norm_mat[px_,px_]=1.0/math_library.sqrt(self.beta_x)

        self.norm_mat[y_,y_]=math_library.sqrt(self.beta_y)
        self.norm_mat[py_,y_]=-self.alpha_y/math_library.sqrt(self.beta_y)
        self.norm_mat[py_,py_]=1.0/math_library.sqrt(self.beta_y)

        self.norm_mat[x_,de_]=self.eta_x
        self.norm_mat[px_,de_]=self.eta_xp
        self.norm_mat[y_,de_]=self.eta_y
        self.norm_mat[py_,de_]=self.eta_yp

    def get_closed_optics_array(self, dpop=0):
        return np.array([self.beta_x.cst()*(1+dpop), self.alpha_x.cst(), self.matched_tune_x.cst(), self.eta_x.cst()*(1+dpop), self.eta_xp.cst(),
                         self.beta_y.cst()*(1+dpop), self.alpha_y.cst(), self.matched_tune_y.cst(), self.eta_y.cst()*(1+dpop), self.eta_yp.cst(),
                         self.c_alpha_c.cst()*(1+dpop), self.c_alpha_c_2.cst()*(1+dpop)*(1+dpop),
                         self.chromaticity_x.cst()*(1+dpop), self.chromaticity_x_2.cst()*(1+dpop)*(1+dpop),
                         self.chromaticity_y.cst()*(1+dpop), self.chromaticity_y_2.cst()*(1+dpop)*(1+dpop),
                         ])

    def get_optics_array(self, dpop=0):
        return np.array([self.beta_x.cst()*(1+dpop), self.alpha_x.cst(), self.tune_x.cst(), self.eta_x.cst()*(1+dpop), self.eta_xp.cst(),
                         self.beta_y.cst()*(1+dpop), self.alpha_y.cst(), self.tune_y.cst(), self.eta_y.cst()*(1+dpop), self.eta_yp.cst()])


def get_matched_optics(one_turn_map, math_library=np, tps6d=None):

    phase_x=math_library.arccos((one_turn_map[x_,x_]+one_turn_map[px_,px_])/2.0)
    beta_x=one_turn_map[x_, px_]/math_library.sin(phase_x)
    if beta_x<0:
        beta_x*=(-1.0)
        phase_x*=2*np.pi-phase_x
    alpha_x=(one_turn_map[x_,x_]-one_turn_map[px_,px_])/2.0/math_library.sin(phase_x)
    eta_x=(one_turn_map[x_, de_]*(1.0-one_turn_map[px_, px_])+one_turn_map[x_, px_]*one_turn_map[px_, de_])/(2.0-one_turn_map[x_,x_]-one_turn_map[px_,px_])
    eta_xp=(one_turn_map[px_, de_]*(1.0-one_turn_map[x_, x_])+one_turn_map[px_, x_]*one_turn_map[x_, de_])/(2.0-one_turn_map[x_,x_]-one_turn_map[px_,px_])

    phase_y=math_library.arccos((one_turn_map[y_,y_]+one_turn_map[py_,py_])/2.0)
    beta_y=one_turn_map[y_, py_]/math_library.sin(phase_y)
    if beta_y<0:
        beta_y*=(-1.0)
        phase_y*=2*np.pi-phase_y
    alpha_y=(one_turn_map[y_,y_]-one_turn_map[py_,py_])/2.0/math_library.sin(phase_y)
    eta_y=(one_turn_map[y_, de_]*(1.0-one_turn_map[py_, py_])+one_turn_map[y_, py_]*one_turn_map[py_, de_])/(2.0-one_turn_map[y_,y_]-one_turn_map[py_,py_])
    eta_yp=(one_turn_map[py_, de_]*(1.0-one_turn_map[y_, y_])+one_turn_map[py_, y_]*one_turn_map[y_, de_])/(2.0-one_turn_map[y_,y_]-one_turn_map[py_,py_])


    temp=Optics(math_library=math_library, beta_x=beta_x, alpha_x=alpha_x, matched_tune_x=phase_x/2/np.pi, eta_x=eta_x, eta_xp=eta_xp,
                beta_y=beta_y, alpha_y=alpha_y, matched_tune_y=phase_y/2/np.pi, eta_y=eta_y, eta_yp=eta_yp)
    if tps6d is not None:
        temp.chromaticity_x=total_derivative(temp.matched_tune_x, 5, tps6d,exclude_list=[4])
        temp.chromaticity_x_2=total_derivative(temp.chromaticity_x, 5, tps6d,exclude_list=[4])
        temp.chromaticity_y=total_derivative(temp.matched_tune_y, 5, tps6d,exclude_list=[4])
        temp.chromaticity_y_2=total_derivative(temp.chromaticity_y, 5, tps6d,exclude_list=[4])
        temp.c_alpha_c=total_derivative(tps6d[4], 5, tps6d, exclude_list=[4])
        temp.c_alpha_c_2=total_derivative(temp.c_alpha_c, 5, tps6d)

    temp.pass_diff=one_turn_map[z_,de_]+temp.H_func_x*math_library.sin(phase_x)*2*(1-math_library.cos(phase_x))+temp.H_func_y*math_library.sin(phase_y)*2*(1-math_library.cos(phase_y))
    return temp


def get_optics(ini_optics, transfer_map, math_library=np):
    #print(transfer_map.dtype, transfer_map.shape, ini_optics.norm_mat.shape, ini_optics.norm_mat.dtype)
    tmap=np.dot(transfer_map, ini_optics.norm_mat)

    dx=tmap[x_, de_]
    dpx=tmap[px_, de_]
    dy=tmap[y_, de_]
    dpy=tmap[py_, de_]

    beta_x=tmap[x_,x_]*tmap[x_,x_]+tmap[x_,px_]*tmap[x_,px_]
    phi_x=math_library.arctan(tmap[x_,px_]/tmap[x_,x_])
    alpha_x=-(math_library.cos(phi_x)*tmap[px_,x_]+math_library.sin(phi_x)*tmap[px_,px_])*math_library.sqrt(beta_x)


    beta_y=tmap[y_,y_]*tmap[y_,y_]+tmap[y_,py_]*tmap[y_,py_]
    phi_y=math_library.arctan(tmap[y_,py_]/tmap[y_,y_])
    alpha_y=-(math_library.cos(phi_y)*tmap[py_,y_]+math_library.sin(phi_y)*tmap[py_,py_])*math_library.sqrt(beta_y)


    temp=Optics(math_library=tpsa, beta_x=beta_x, alpha_x=alpha_x, tune_x=phi_x/2/np.pi+ini_optics.tune_x, eta_x=dx, eta_xp=dpx,
                beta_y=beta_y, alpha_y=alpha_y, tune_y=phi_y/2/np.pi+ini_optics.tune_y, eta_y=dy, eta_yp=dpy)
    return temp



