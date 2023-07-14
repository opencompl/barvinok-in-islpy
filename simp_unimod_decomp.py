import numpy as np
from numpy.linalg import inv, det
from islpy import BasicSet, Space, Constraint, Context
from fpylll import LLL, IntegerMatrix
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
        """
        r = np.array(self.rays)
        R = (inv(r) * det(r)).round().astype(int).tolist()
        B = olll.reduction(R, 0.75)
        U = (np.array(B) @ r).tolist()
        #assert ((np.array(U) @ np.array(self.rays) == np.array(B)).all())
        U = [[int(i/gcd(*l)) for i in l] for l in U]

        index, λ = min(enumerate(B), key = lambda r : abs(max(r[1], key=abs)))
        v = U[index]

        if all(map(lambda x : x <= 0, λ)):
            v = [-i for i in v]
        
        return v

    def __repr__(self):
        return f"{'+' if self.sign == 1 else '-'}{self.rays}"

def unimodular_decomp(cone):
    ind = cone.get_index()
    if ind == 1: return [cone]
    else:
        cones = []
        rays = cone.rays
        w = cone.get_sample_point()
        for i in range(len(rays)):
            replaced = subAtWith(rays, i, w)
            d = det(replaced).item()
            if d == 0: continue
            ki = Cone(replaced, sign(d) * cone.sign)
            cones.append(ki)
            print(f"{cone} -> {ki}")
        final_nested = map(unimodular_decomp, cones)
        final = [cone for decomp in final_nested for cone in decomp]
        return final