import numpy as np
from islpy import BasicSet, Space, Constraint, Context
import islpy as isl

# For simplicial cones
class Cone():
    def __init__(self, vertex, rays, sign=1):
        #              [d]     [d, d]
        self.vertex = vertex
        self.rays = rays
        self.sign = sign

    def get_index(self):
        """Volume of parallelepiped"""
        det = np.linalg.det(self.rays)
        return (0.5 * det)

    def get_sample_point(self): # TODO
        """
        Use LLL and find
        """
        pass

def unimodular_decomp(cone):
    ind = cone.get_index()
    if ind == 1: return [(1, cone)]
    else:
        cones = []
        rays = cone.rays
        w = cone.get_sample_point()
        sign = 1                     # TODO
        for i in range(len(rays)):
            ki = Cone(cone.vertex, rays[:i] + [w] + rays[i+1:], sign)
            cones.append(ki)
        final_nested = map(unimodular_decomp, cones)
        final = [cone for decomp in final_nested for cone in decomp]
        return final