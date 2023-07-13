import numpy as np
from islpy import BasicSet, Space, Constraint, Context
from fpylll import LLL, IntegerMatrix

# For simplicial cones
class Cone():
    def __init__(self, rays, sign=1):
        #              [d, d]
        self.rays = rays
        self.sign = sign
        self.d = len(self.rays)

    def get_index(self):
        """Volume of parallelepiped"""
        det = abs(np.linalg.det(self.rays))
        return (0.5 * det)

    def get_sample_point(self):
        """
        Use LLL and find
        """
        U = IntegerMatrix.identity(self.d)
        B = IntegerMatrix.from_matrix(self.rays)
        
        LLL.reduction(B, U)
        B = [list(r) for r in list(B)]
        U = [list(r) for r in list(U)]

        index, λ = min(enumerate(B), key = lambda r : max(r[1]))
        v = U[index]

        if all(map(lambda x : x <= 0, λ)):
            λ = -λ
            v = -v
        
        firstNonzeroElem = [x for x in λ if x != 0][0]
        sign = firstNonzeroElem/abs(firstNonzeroElem)
        
        return v, sign

def unimodular_decomp(cone):
    ind = cone.get_index()
    if ind == 1: return [cone]
    else:
        cones = []
        rays = cone.rays
        w, sign = cone.get_sample_point()
        for i in range(len(rays)):
            ki = Cone(rays[:i] + [w] + rays[i+1:], sign * cone.sign)
            cones.append(ki)
        final_nested = map(unimodular_decomp, cones)
        final = [cone for decomp in final_nested for cone in decomp]
        return final