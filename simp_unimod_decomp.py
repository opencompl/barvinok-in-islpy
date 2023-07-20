import numpy as np
from numpy.linalg import inv, det
import olll
from math import gcd
import cdd

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
        d = abs(det(self.rays).round().astype(int))
        return d

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

def triangulate(cone):
    """
    Use Delaunay's method:
    * add a coefficient which is the sum of squares of the other coefficients
    * find the facets of the polyhedron formed by these points
    * filter them according as the last coordinate in their outer normal vector is negative
    * project these facets down
    """
    extended_rays = []
    for r in cone.rays:
        extended_rays.append(r + [sum(x**2 for x in r)])

    mat = cdd.Matrix(extended_rays)
    mat.rep_type = cdd.RepType.GENERATOR
    lifted_poly = cdd.Polyhedron(mat)
    # This polyhedron is in V (generator) representation

    facets_of_lifted_poly = lifted_poly.get_inequalities()
    lower_facet_indices = []
    for i in range(facets_of_lifted_poly.col_size):
        if facets_of_lifted_poly[i][-1] < 0:
            lower_facet_indices.append(i) # last coordinate of facet
    
    facets_to_incident_vertices = lifted_poly.get_incidence()
    # Since in V-representation, we get F-to-V map

    triangle_generator_indices = [facets_to_incident_vertices[i] for i in lower_facet_indices]
    # (Indices of) generators of lower facets

    triangles = []
    for generator_index_set in triangle_generator_indices:
        generator_index_list = list(generator_index_set)
        triangles.append(Cone([cone.rays[i] for i in generator_index_list], cone.sign))
    
    return triangles

"""Tests
>>> unimodular_decomp(Cone([[-3,1,1],[-1,-3,-1],[-1,-2,-1]]))
[+[[-1, 0, 0], [-1, -3, -1], [-1, -2, -1]], -[[-3, 1, 1], [-1, 0, 0], [-1, -2, -1]], +[[-2, -1, 0], [-1, -3, -1], [-1, 0, 0]], +[[-3, 1, 1], [-2, -1, 0], [-1, 0, 0]]]
>>> unimodular_decomp(Cone([[-2,-3,0], [-1,0,-2], [-1,0,0]]))
[+[[-1, -1, 0], [0, 0, -1], [-1, 0, 0]], -[[-1, -1, 0], [-1, 0, -2], [0, 0, -1]], +[[-2, -3, 0], [-1, 0, -2], [-1, -1, 0]]]

"""