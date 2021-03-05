from .. element import Element
from .. track import pass_components
import copy
class Drift(Element):
    propertyNames = copy.deepcopy(Element.propertyNames)
    propertyNames.update(['EXACT_DRIFT'])
    def __init__(self, name, **param):
        Element.__init__(self, name, self.__class__.__name__)
        self.properties['EXACT_DRIFT'] = True
        self.set_param(**param)

    def track(self, coor, math_lib):
        pass_components.pass_drift(coor, self.properties['L'], math_lib=math_lib)

        # def build_CL_program(self, ctx):
        #    Drift.pass_cl_kernels= cl.Program(ctx, CL_track.kernel_pass_drift).build()

        # def track_GPU(self, coor_buf, ctx, queue, num_p, dim):
        #    Drift.pass_cl_kernels.pass_drift.set_scalar_arg_dtypes([None, np.float, np.uint32,np.uint32])
        #    Drift.pass_cl_kernels.pass_drift(coor_buf, self.properties['L'], num_p, dim)
