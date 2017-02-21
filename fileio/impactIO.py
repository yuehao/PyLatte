from .fileIO import LatticeFile
import numpy as np

class ImpactInternalDistribution(object):
    def __init__(self, ref_freq=100.0e6, Ek0=1.0e6, mass=510.99e3, npart=[1,], charge=[1.0,], charge_state=[1.0/510.99e3,]):
        self.vgamma=(1+Ek0)/mass
        self.vbeta=np.sqrt(1.0-1.0/self.vgamma/self.vgamma)
        self.p0c=self.vgamma*self.vbeta*mass

        self.set_charge(npart, charge, charge_state)

    def set_charge(self, npart=[1,], charge=[1.0,], charge_state=[1.0/511.0e3,]):

        if isinstance(npart, int) or isinstance(npart, float):
            self.npart=[int(npart), ]
            self.total_npart=np.sum(self.npart)
            if isinstance(charge, float) or isinstance(charge, int):
                self.charge = np.array([charge, ])
            else:
                self.charge = np.array([charge[0], ])

            if isinstance(charge_state, float):
                self.charge_state = np.array([charge_state, ])
            else:
                self.charge_state = np.array([charge_state[0], ])



        elif isinstance(npart, list) or isinstance(npart, np.ndarray):
            self.npart = np.array(npart)
            self.total_npart = np.sum(self.npart)

            if len(self.npart) <= len(charge_state):
                self.charge_state = np.array(charge_state)[0:len(self.npart)]
            else:
                print('Data for charge_state is not enough')

            if len(self.npart) <= len(charge):
                self.charge = np.array(charge)[0:len(self.npart)]
            else:
                print('Data for charge is not enough')

        self.coordinate = np.zeros((self.total_npart, 9))
        self.coordinate[:, -1] = np.arange(int(self.total_npart)) + 1
        if len(self.npart) == 1:
            self.coordinate[:, -2] += self.charge[0]
            self.coordinate[:, -3] += self.charge_state[0]
        else:
            split = np.cumsum(self.npart)[:-1]
            temp = np.split(self.coordinate[:, -2])
            temp += self.charge
            self.coordinate[:, -2] = np.concatenate(temp)

            temp = np.split(self.coordinate[:, -3])
            temp += self.charge
            self.coordinate[:, -3] = np.concatenate(temp)

    def write(self, filename='partcl.data'):
        with open(filename, 'w') as f:
            f.write('{}\t0\t0\n'.format(int(self.total_npart)))
            for i in range(int(self.total_npart)):
                for j in range(9):
                    f.write('{}\t'.format(self.coordinate[i,j]))
                f.write('\n')

    def eleven_points(self, d=1e-4,dp=1e-6,de=0, offset=np.zeros(6), charge=-1.0, charge_state=1.0/511e3):
        self.set_charge(npart=11, charge=charge, charge_state=charge_state)
        self.coordinate[:,:6] = np.array(offset)
        self.coordinate[1, 0] += d
        self.coordinate[2, 0] += -d
        self.coordinate[3, 1] += dp
        self.coordinate[4, 1] += -dp
        self.coordinate[5, 2] += d
        self.coordinate[6, 2] += -d
        self.coordinate[7, 3] += dp
        self.coordinate[8, 3] += -dp
        self.coordinate[9, 5] += de
        self.coordinate[10, 5] += -de




class ImpactHeader(object):
    def __init__(self):
        # line 1
        self.external_distribution=None
        self.ncpu = 1

        # line 2
        self.n_dim = 6
        self.n_part = 1
        self.integrator = 2  # 1 for linear; 2 for nonlinear
        self.error_study = False
        self.output_type = 1  # 1 for std; 2 for %95

        # line_3
        self.grid = [64, 64, 64]
        self.boundary_type = 1  # 1 open 3D; 2 open xy, z peroidic; 3 round, z open; 4 round, z periodic; 5 rect, z open; 6 rect, z periodic
        self.grid_boundary = [0.01, 0.01, 0.1]

        # line_4
        self.distribution = 3
        # 1 Rect; 2 Gaussian; 3 Waterbag; 4 Uniform X, Gaussian P; 5 KV Tran Uniform Long;
        # 16 MC Waterbag; 17 MC Gaussian; 19 read from 'partcl.data'
        self.restart = 0
        self.sub_cycle = 0
        self.n_charge_state = 1

        # line_5
        self.n_part_mc = [1, ]

        # line_6
        self.current_mc = [0.0, ]

        # line_7, Q/m 1/(eV)
        self.charge_over_mass = [1.0 / 9.3929e8, ]

        # line_8-10 sigmax    lambdax      mux     mismatchx  mismatchpx offsetX  offsetPx
        self.moment_x = [0.22734189E-02, 0.88312578E-04, 0.00000000E+00, 1.000, 1.000, 0.000, 0.000]
        self.moment_y = [0.22734189E-02, 0.88312578E-04, 0.00000000E+00, 1.000, 1.000, 0.000, 0.000]
        self.moment_z = [0.76704772E-01, 0.34741445E-05, 0.00000000E+00, 1.000, 1.000, 0.000, 0.000]

        # line_11 current_ave, ini_kinetic_energy, mass_in_eV, charge,freq,iniphase
        self.ave_current = 0.001
        self.ini_kenergy = 100e3
        self.mass = 9.3929e8
        self.charge = 1.0
        self.freq = 80.5e6
        self.ini_phase = 0.0

    def check_boundary(self):
        import math
        def grid_open(i):
            return int(math.pow(2, math.floor(math.log(i, 2))))

        if self.boundary_type == 1:
            self.grid[0] = grid_open(self.grid[0])
            self.grid[1] = grid_open(self.grid[1])
            self.grid[2] = grid_open(self.grid[2])
        elif self.boundary_type == 2:
            self.grid[0] = grid_open(self.grid[0])
            self.grid[1] = grid_open(self.grid[1])
            self.grid[2] = grid_open(self.grid[2]) + 1
        elif self.boundary_type == 3:
            self.grid[0] = grid_open(self.grid[0])
            self.grid[1] = grid_open(self.grid[1]) + 1
            self.grid[2] = grid_open(self.grid[2])
        elif self.boundary_type == 4:
            self.grid[0] = grid_open(self.grid[0])
            self.grid[1] = grid_open(self.grid[1]) + 1
            self.grid[2] = grid_open(self.grid[2]) + 1
        elif self.boundary_type == 5:
            self.grid[0] = grid_open(self.grid[0]) + 1
            self.grid[1] = grid_open(self.grid[1]) + 1
            self.grid[2] = grid_open(self.grid[2])
        elif self.boundary_type == 6:
            self.grid[0] = grid_open(self.grid[0]) + 1
            self.grid[1] = grid_open(self.grid[1]) + 1
            self.grid[2] = grid_open(self.grid[2]) + 1
        else:
            print("Undefined boundary condition")
            exit(-1)
    def check(self):
        if self.distribution == 19:
            assert isinstance(self.external_distribution, ImpactInternalDistribution)
            self.external_distribution.write()
            self.n_part=self.external_distribution.total_npart

        assert len(self.n_part_mc) == self.n_charge_state
        assert len(self.charge_over_mass) == self.n_charge_state
        assert len(self.current_mc) == self.n_charge_state
        assert len(self.moment_x) == 7
        assert len(self.moment_y) == 7
        assert len(self.moment_z) == 7
        if self.n_charge_state == 1:
            self.charge_over_mass[0]=np.abs(self.charge/self.mass)


    def set_param(self, **param_dict):
        namelist = self.__dict__.keys()
        for k, v in param_dict.items():
            if k in namelist:
                self.__dict__[k] = v
                if k == 'boundary_type':
                    self.check_boundary()
            else:
                print("Warning, {} is not a recognized key, ignored the action {}={}.".format(k, k, v))
                exit(-1)
        self.check()

    def write(self, f):
        #line1
        f.write('{ncpu1}\t{ncpu2}\n'.format(ncpu1=self.ncpu, ncpu2=1))
        #line2
        f.write('{nd}\t{nmp}\t{intg}\t{err}\t{out}\n'.format(nd=self.n_dim, nmp=self.n_part, intg=self.integrator,
                                                           err=int(self.error_study), out=self.output_type))
        #line3
        f.write('{g[0]}\t{g[1]}\t{g[2]}\t{bound}\t{be[0]}\t{be[1]}\t{be[2]}\n'.format(g=self.grid, bound=self.boundary_type,
                                                                                    be=self.grid_boundary))
        #line4
        f.write('{dis}\t{rst}\t{sc}\t{cs}\n'.format(dis=self.distribution, rst=self.restart,
                                                  sc=self.sub_cycle, cs=self.n_charge_state))
        #line5
        for i in range(self.n_charge_state):
            f.write('{}\t'.format(self.n_part_mc[i]))
        f.write('\n')

        # line6
        for i in range(self.n_charge_state):
            f.write('{}\t'.format(self.current_mc[i]))
        f.write('\n')

        # line7
        for i in range(self.n_charge_state):
            f.write('{}\t'.format(self.charge_over_mass[i]))
        f.write('\n')

        # line8-10
        for temp in (self.moment_x):
            f.write('{}\t'.format(temp))
        f.write('\n')
        for temp in (self.moment_y):
            f.write('{}\t'.format(temp))
        f.write('\n')
        for temp in (self.moment_z):
            f.write('{}\t'.format(temp))
        f.write('\n')

        #line11
        f.write('{cur}\t{ke}\t{mass}\t{ch}\t{f}\t{ph}\t99.9\n'.format(cur=self.ave_current, ke=self.ini_kenergy, mass=self.mass,
                                                                ch=self.charge, f=self.freq, ph=self.ini_phase))


class ImpactInputFile(LatticeFile):
    element_type={'DRIFT':0,
                  'QUADRUPOLE':1,
                  'CONSTFOCUS':2,
                  'SOLENOID':3,
                  'DIPOLE':4,
                  'DTL':101,
                  'CCDTL':102,
                  'CCL':103,
                  'CAVITY':104,
                  'SOLENOIDRF':105,
                  'EMFIELD':110,
                  'CENTER': -1,
                  'BPM': -2,
                  'PROFILE': -2,
                  }
    common_parameter={'L':0.0, 'N_STEP':10, 'MAP_STEP':20, 'APERTURE':0.01, 'FILEID':0}
    error_parameter={'DX':0.0, 'DY':0.0, 'DANGX':0.0, 'DANGY':0.0, 'DANGZ':0}
    parameters={
                'QUADRUPOLE': {'GRADIENT':0.0, },
                'CONSTFOCUS': {'KX':0.0, 'KY':0.0, 'KZ':0.0,} ,
                'SOLENOID': {'BZ':0.0,},
                'DIPOLE':{'ANGLE':0.0,'K1':0.0,'E1':0.0, 'E2':0.0, 'H1':0.0, 'H2':0.0, 'FINT':0.0},
                'DTL': {'SCALE':0.0, 'FREQ':0.0, 'PHASE':0.0, 'Q1L':0.0, 'Q1G':0.0, 'Q2L':0.0, 'Q2G':0.0},
                'CCDTL': {'SCALE':0.0, 'FREQ':0 ,'PHASE':0.0},
                'CCL': {'SCALE':0.0, 'FREQ':0 ,'PHASE':0.0},
                'CAVITY': {'SCALE':0.0, 'FREQ':0 ,'PHASE':0.0},
                'SOLENOIDRF': {'SCALE':0.0, 'FREQ':0.0,'PHASE':0.0, 'BZ':0.0},
                'EMFIELD': {'SCALE':0.0, 'FREQ':0.0,'PHASE':0.0, 'BZ':0.0, 'MODE':1.0, 'COOR':1.0},

                }
    def __init__(self):
        LatticeFile.__init__(self)
        self.header = ImpactHeader()

    def write(self, filename='test.in'):
        with open(filename, 'w') as f:
            self.header.write(f)
            if self.useline == '':
                print('No line is defined. need to write manually')
                return
            ind=self.beamlineNameDict[self.useline]

            for elename in self.beamlineList[ind]['LINE']:
                ele=self.elementList[self.getElementIndex(elename)]
                dtemp={}
                dtemp=ImpactInputFile.common_parameter+ImpactInputFile.error_parameter+ImpactInputFile.parameters[ele['TYPE']]
                dtemp.update(ele)
                if ele['TYPE']=='DRIFT':
                    f.write('{l}\t{s1}\t{s2}\t{tp}\t{r}\t/'.format(l=dtemp['L'],s1=dtemp['STEPS'], s2=dtemp['MAPSTEPS'],
                                                                tp=ImpactInputFile.element_type[ele['TYPE']],r=dtemp['RADIUS']))
                if ele['TYPE']=='QUADRUPOLE':
                    f.write('{l}\t{s1}\t{s2}\t{tp}\t{g}\t{id}\t{r}\t{dx}\t{dy}\t{ax}\t{ay}\t{az}/'.
                            format(l=dtemp['L'],s1=dtemp['N_STEP'], s2=dtemp['MAP_STEP'], tp=ImpactInputFile.element_type[ele['TYPE']],
                                   r=dtemp['RADIUS'],g=dtemp['GRADIENT'],id=dtemp['FILEID'], dx=dtemp['DX'],dy=dtemp['DY'],
                                   ax=dtemp['DANGX'],ay=dtemp['DANGY'],az=dtemp['DANGZ']))
                if ele['TYPE'] == 'CONSTFOCUS':
                    f.write('{l}\t{s1}\t{s2}\t{tp}\t{kx}\t{ky}\t{kz}\t{r}/'.
                            format(l=dtemp['L'], s1=dtemp['STEPS'], s2=dtemp['MAPSTEPS'],tp=ImpactInputFile.element_type[ele['TYPE']],
                                   r=dtemp['RADIUS'], kx=dtemp['KX'], ky=dtemp['KY'], kz=dtemp['KZ']))
                if ele['TYPE']=='SOLENOID':
                    f.write('{l}\t{s1}\t{s2}\t{tp}\t{bz0}\t{id}\t{r}\t{dx}\t{dy}\t{ax}\t{ay}\t{az}/'.
                            format(l=dtemp['L'],s1=dtemp['STEPS'], s2=dtemp['MAPSTEPS'], tp=ImpactInputFile.element_type[ele['TYPE']],
                                   r=dtemp['RADIUS'],bz0=dtemp['GRADIENT'],id=dtemp['FILEID'], dx=dtemp['DX'],dy=dtemp['DY'],
                                   ax=dtemp['DANGX'],ay=dtemp['DANGY'],az=dtemp['DANGZ']))





    def checkType(self, typename, parameterName=None):
        if typename.upper() in ImpactInputFile.element_type:
            if parameterName is None:
                return True
            else:
                parameter_dict={}
                parameter_dict.update(ImpactInputFile.common_parameter)
                parameter_dict.update(ImpactInputFile.error_parameter)
                parameter_dict.update(ImpactInputFile.parameters[typename.upper()])
                return parameterName in parameter_dict





class ImpactOutput(object):
    '''
    Regroup the structure to
    self.s , self.energy, self.vbeta self.aperture (from fort.18)
    self.moment_x/y/z: 0:cx,1:cxp,2:sx,3:sxp,4:norm.emit,5:max_x,6:max_px,7:m3_x,8:m3_px,9:m4_x,10:m4_px,
    self.twiss: 0:betax, 1:alphax, 2:betay, 3:alphay
    '''
    def __init__(self):
        import glob
        self.found_files=glob.glob("fort.*")
        if not self.found_files:
            print('Cannot find any outfile in this directory, Exiting\n')
            exit(-1)

    def read_files(self):
        for fn in self.found_files:
            otype=int(fn.split('.')[-1])
            if otype == 18:
                temp=np.genfromtxt(fn)
                self.s = temp[:, 0]
                self.vgamma = temp[:, 2]
                self.energy = temp[:, 3]
                self.vbeta = temp[:, 4]
                self.aperture = temp[:, 5]
                self.moment_x = np.zeros((len(self.s), 11))
                self.moment_y = np.zeros((len(self.s), 11))
                self.moment_z = np.zeros((len(self.s), 11))

                self.twiss = np.zeros((len(self.s), 4))

            elif otype == 24:
                temp = np.genfromtxt(fn)
                self.moment_x[:, 0] = temp[:, 1]  # ave_x
                self.moment_x[:, 1] = temp[:, 3]  # ave_xp
                self.moment_x[:, 2] = temp[:, 2]  # sig_x
                self.moment_x[:, 3] = temp[:, 4]  # sig_xp
                self.moment_x[:, 4] = temp[:, 6]  # emittance

                self.twiss[:, 0]=temp[:, 6]
                self.twiss[:, 1]=temp[:, 5]

            elif otype == 25:
                temp = np.genfromtxt(fn)
                self.moment_y[:, 0] = temp[:, 1]
                self.moment_y[:, 1] = temp[:, 3]
                self.moment_y[:, 2] = temp[:, 2]
                self.moment_y[:, 3] = temp[:, 4]
                self.moment_y[:, 4] = temp[:, 6]

                self.twiss[:, 2] = self.moment_y[:, 2] * self.moment_y[:, 2] / self.moment_y[:, 4]
                self.twiss[:, 3] = temp[:, 5]

            elif otype == 26:
                temp = np.genfromtxt(fn)
                self.moment_z[:, 0] = temp[:, 1]
                self.moment_z[:, 1] = temp[:, 3]
                self.moment_z[:, 2] = temp[:, 2]
                self.moment_z[:, 3] = temp[:, 4]
                self.moment_z[:, 4] = temp[:, 6]

            elif otype == 27:
                temp = np.genfromtxt(fn)
                self.moment_x[:, 5:7] = temp[:, 1:3]
                self.moment_y[:, 5:7] = temp[:, 3:5]

            elif otype == 29:
                temp = np.genfromtxt(fn)
                self.moment_x[:, 7:9] = temp[:, 1:3]
                self.moment_y[:, 7:9] = temp[:, 3:5]

            elif otype == 30:
                temp = np.genfromtxt(fn)
                self.moment_x[:, 9:] = temp[:, 1:3]
                self.moment_y[:, 9:] = temp[:, 3:5]

            elif otype == 32:
                self.nparticle = np.genfromtxt(fn)[:,1]








