import numpy as np
class field3D(object):
    def __init__(self, filename, header=3):
        self.read_file(filename, header)
    def read_file(self, filename, header, header_seq=['x','y','z']):
        with open(filename, 'r') as f:
            for dim in header_seq:
                l=f.readline().split()
                self.__dict__['{}range'.format(dim)]=np.array(l)[0:2].astype(float)
                self.__dict__['{}size'.format(dim)]=int(l[-1])+1
                self.__dict__['{}coor'.format(dim)] = np.linspace(float(l[0]),float(l[1]),int(l[2])+1)
        self.field=np.genfromtxt(filename, skip_header=header)



    def interpolate_all(self, zbin, xc, yc, nint=10):

        slice_num = self.xsize * self.ysize
        xmin = self.xsize / 2 - nint / 2;
        xmax = self.xsize / 2 + nint / 2 + 1;
        ymin = self.ysize / 2 - nint / 2;
        ymax = self.ysize / 2 + nint / 2 + 1;

        xcoefd = np.zeros_like(self.xcoor[xmin:xmax])
        ycoefd = np.zeros_like(self.ycoor[ymin:ymax])

        for i in range(len(xcoefd)):
            tempx = np.roll(self.xcoor[xmin:xmax], -i)
            xcoefd[i] = np.prod((tempx[0] - tempx)[1:])
        for i in range(len(ycoefd)):
            tempy = np.roll(self.ycoor[ymin:ymax], -i)
            ycoefd[i] = np.prod((tempy[0] - tempy)[1:])
        slice = self.field[slice_num * zbin:slice_num * (zbin + 1), :].reshape((self.ysize,self.xsize,6))[ymin:ymax,xmin:xmax,:]

        xcoefn = np.zeros_like(self.xcoor[xmin:xmax])+xc
        ycoefn = np.zeros_like(self.ycoor[ymin:ymax])+yc

        for i in range(len(xcoefn)):
            tempx = np.roll(self.xcoor[xmin:xmax], -i)
            terms= -(tempx[1:])+xc
            xcoefn[i] = np.prod(terms)

        for i in range(len(ycoefn)):
            tempy = np.roll(self.ycoor[ymin:ymax], -i)
            terms = - tempy[1:] + yc
            ycoefn[i] = np.prod(terms)

        return xcoefn












