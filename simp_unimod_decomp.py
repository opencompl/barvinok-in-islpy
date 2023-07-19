import numpy as np
from numpy.linalg import inv, det
from islpy import BasicSet, Space, Constraint, Context
import olll
from math import gcd

def sign(x): return x/abs(x)

def subAtWith(l, i, x): return (l[:i] + [x] + l[i+1:])

# For simplicial cones
class Cone():
    def __init__(self, rays, sign=1):
        #              [d, d]
        self.rays = rays
        self.sign = sign
        self.d = len(self.rays)

    def get_index(self):
        """Volume of parallelepiped"""
        det = int(abs(np.linalg.det(self.rays)))
        return det

    def get_sample_point(self):
        """
        Use LLL and find
        https://math.ucdavis.edu/~deloera/researchsummary/barvinokalgorithm-latte1.pdf, p. 1279-80
        """
        r = np.array(self.rays)
        R = (inv(r) * det(r)).round().astype(int).tolist()
        B = olll.reduction(R, 0.75)
        U = (np.array(B) @ r).tolist()
        #assert ((np.array(U) @ np.array(self.rays) == np.array(B)).all())
        U = [[int(i/gcd(*l)) for i in l] for l in U]

        index, λ = min(enumerate(B), key = lambda r : max(map(abs,r[1])))
        v = U[index]

        if all(map(lambda x : x <= 0, λ)):
            v = [-i for i in v]
            λ = [-i for i in λ]
        
        return v, λ

    def __repr__(self):
        return f"{'+' if self.sign == 1 else '-'}{self.rays}"

def unimodular_decomp(cone):
    ind = cone.get_index()
    assert(ind != 0)
    if ind == 1: return [cone]
    else:
        cones = []
        rays = cone.rays
        w, λ = cone.get_sample_point()
        for i in range(len(rays)):
            if λ[i] == 0: continue
            replaced = subAtWith(rays, i, w)
            ki = Cone(replaced, sign(λ[i]) * cone.sign)
            cones.append(ki)
        final_nested = map(unimodular_decomp, cones)
        final = [cone for decomp in final_nested for cone in decomp]
        return final
