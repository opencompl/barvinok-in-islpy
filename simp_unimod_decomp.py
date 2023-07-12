import numpy as np
from islpy import BasicSet, Space, Constraint, Context
import islpy as isl

# For simplicial cones
class Cone():
    def __init__(self, coeffs, constants):
        #              [d, d]  [d]
        self.coeffs = coeffs
        self.constants = constants
        self.d = len(coeffs)

        space = Space.set_alloc(Context(), nparam=0, dim=self.d)
        set = BasicSet.universe(space)

        for i in range(self.d):
            cons = Constraint.alloc_inequality(space)
            cons = cons.set_constant_val(self.constants[i])
            for j, c in enumerate(self.coeffs[i]):
                cons = cons.set_coefficient_val(isl.dim_type.set, j, c)
                set = set.add_constraint(cons)

        self.set = set

    def get_index(self):
        """Volume of parallelepiped"""
        det = np.linalg.det(self.rays)
        return (0.5 * det)

    def get_rays(self): # TODO
        pass

def unimodular_decomp(cone):
    ind = cone.get_index()
    if ind == 1: return [cone]
    else:
        cones = []
        rays = cone.get_rays()
        w = cone.set.sample_point()
        for i in range(len(rays)):
            ki = Cone(rays[:i] + [w] + rays[i+1:])   # TODO
            cones.append(ki)
        final_nested = map(unimodular_decomp, cones)
        final = [cone for decomp in final_nested for cone in decomp]
        return final