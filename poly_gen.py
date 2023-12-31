import numpy as np
from numpy.linalg import inv, det
from numpy import transpose, ceil
import olll
from math import gcd, factorial
import cdd
from functools import reduce

def sign(x): return x/abs(x)

def subAtWith(l, i, x): return (l[:i] + [x] + l[i+1:])

def dot(a, b): return sum(ai*bi for ai,bi in zip(a,b))

def get_index(cone):
    """Volume of parallelepiped"""
    rays = [row[1:] for row in cone.get_generators()]
    if len(rays) > len(rays[0]): return 0
    # There are more rays than dimensions of the space
    d = abs(det(rays).round().astype(int))
    return d

def get_sample_point(cone):
    """
    Use LLL and find
    https://math.ucdavis.edu/~deloera/researchsummary/barvinokalgorithm-latte1.pdf, p. 1279-80
    """
    rays = [row[1:] for row in cone.generators()]
    r = np.array(rays)
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

def get_dual(cone):
    mat = cone.get_generators()
    mat.rep_type = cdd.RepType.INEQUALITY
    dual = cdd.Polyhedron(mat)
    return dual

def unimodular_decomp(cone):
    dual = get_dual(cone)
    ind = get_index(dual)
    if ind == 0:
        simplicial = triangulate(dual)
        dual_final_nested = [unimodular_decomp_simplicial(1, c) for c in simplicial]
        dual_final = [c for l in dual_final_nested for c in l]
    else:
        dual_final = unimodular_decomp_simplicial(1, dual)
    final = [(s, get_dual(c)) for s, c in dual_final]
    return final

def unimodular_decomp_simplicial(s, cone):
    ind = get_index(cone)
    if ind == 1: return [(1, cone)]
    else:
        cones = []
        rays = cone.get_generators()
        w, λ = get_sample_point(cone)
        for i in range(len(rays)):
            if λ[i] == 0: continue
            replaced = subAtWith(rays, i, [0] + w)
            mat = cdd.Matrix(replaced)
            mat.rep_type = cdd.RepType.GENERATOR
            ki = cdd.Polyhedron(mat)
            cones.append((sign(λ[i]) * s, ki))
        final_nested = [unimodular_decomp_simplicial(s, c) for s, c in cones]
        final = [s_cone for decomp in final_nested for s_cone in decomp]
        return final

def triangulate(cone):
    """
    Use Delaunay's method:
    * add a coefficient which is the sum of squares of the other coefficients
    * find the facets of the polyhedron formed by these points
    * filter them according as the last coordinate in their outer normal vector is negative ("lower facets")
    * project these facets down
    """
    rays = cone.get_generators()
    extended_rays = []
    for r in rays:
        extended_rays.append(r + [sum(x**2 for x in r)])

    mat = cdd.Matrix(extended_rays)
    mat.rep_type = cdd.RepType.GENERATOR
    lifted_poly = cdd.Polyhedron(mat)
    # This polyhedron is in V (generator) representation

    facets_of_lifted_poly = lifted_poly.get_inequalities()
    lower_facet_indices = []
    for i in range(facets_of_lifted_poly.row_size):
        if facets_of_lifted_poly[i][-1] < 0:
            lower_facet_indices.append(i) # last coordinate of facet
    
    facets_to_incident_vertices = lifted_poly.get_incidence()
    # Since in V-representation, we get F-to-V map

    triangle_generator_indices = [facets_to_incident_vertices[i] for i in lower_facet_indices]
    # (Indices of) generators of lower facets

    triangles = []
    for generator_index_set in triangle_generator_indices:
        generator_index_list = list(generator_index_set)
        mat = cdd.Matrix([rays[i] for i in generator_index_list])
        mat.rep_type = cdd.RepType.GENERATOR
        triangles.append(cdd.Polyhedron(mat))
        # Get the lower-dimensional projections of the generators for each facet
    
    return triangles

class GenFunc():
    def __init__(self, den_exps : list[list[list[int]]], num_exps : list[list[int]], signs : list[int]):
        self.den_exps = den_exps
        self.num_exps = num_exps
        self.signs = signs
    
    def add(self, other): # concat the lists, one element per term
        self.den_exps += other.den_exps
        self.num_exps += other.num_exps
        self.signs += other.signs

    def __getitem__(self, i):
        return self.signs[i], self.num_exps[i], self.den_exps[i]
    
    def __len__(self):
        return len(self.signs)
    
    def __repr__(self):
        return "\n".join([
                    f"  x^{list(n)}\n" +
                    f"{'+' if s == 1 else '-'} ---------\n  " +
                    "".join([f"(1 - x^{list(d)})" for d in ds]) + "\n"
                        for s, n, ds in zip(self.signs, self.num_exps, self.den_exps)])

def unimod_cone_gen_func(vertex, s_cone):
    s, cone = s_cone
    if all(isinstance(i, int) for i in vertex): numerator = vertex
    else:
        λ = inv(transpose(cone.rays)) @ vertex
        numerator = ceil(λ).astype(int)
    
    denominator = [row[1:] for row in cone.get_generators()]

    return GenFunc(den_exps = [denominator],
                   num_exps = [numerator],
                   signs = [s])

def polytope_gen_func(poly): # polyhedron in H-rep (Ax <= b)
    vertices = poly.get_generators()
    faces = poly.get_inequalities()
    vertices_to_incident_faces = poly.get_incidence()

    gf = GenFunc([], [], [])
    for i in range(vertices.row_size):
        vertex = vertices[i][1:]

        incident_faces = vertices_to_incident_faces[i]
        cone_matrix_rows = []
        for f in incident_faces:
            cone_matrix_rows.append(faces[f])
        cone_matrix = cdd.Matrix(cone_matrix_rows)
        cone_matrix.rep_type = cdd.RepType.INEQUALITY
        cone = cdd.Polyhedron(cone_matrix)
        rays_and_vertex = cone.get_generators()
        
        rays = [row for row in rays_and_vertex if row[0] == 0]
        mat = cdd.Matrix(rays)
        mat.rep_type = cdd.RepType.GENERATOR

        C = cdd.Polyhedron(mat)

        decomp = unimodular_decomp(C)

        cone_gf = GenFunc([],[],[])
        for unimod_cone in decomp:
            cone_gf.add(unimod_cone_gen_func(vertex, unimod_cone))

        gf.add(cone_gf)
    
    return gf

def comb(n, r):
    numerator = reduce(lambda a, b: a*b, [n-i for i in range(0, r)], 1)
    denominator = factorial(r)
    return numerator/denominator

def conv(l1, l2):
    if len(l1) > len(l2): l2 += [0 for _ in range(len(l1)-len(l2))]
    elif len(l2) > len(l1): l1 += [0 for _ in range(len(l2)-len(l1))]

    l = len(l1)

    return [sum(l1[i] * l2[j-i] for i in range(0,j+1))
                for j in range(0,l)]

def get_coeffs_numerator(num):
    def p_of_s(i):
        return comb(num, i)
    return p_of_s

def get_coeffs_denominator(dens):
    coeffs_per_den = [[comb(den+1, k+1) for k in range(0, den+1)] for den in dens]
    coeffs = reduce(conv, coeffs_per_den[1:], coeffs_per_den[0])
    def q_of_s(i):
        try: return coeffs[i]
        except IndexError: return 0
    return q_of_s

def get_coeff_of(r, p, q):
    if r == 0: return p(0)/q(0)
    else:
        t = p(r)
        for i in range(1,r+1):
            t -= q(i) * get_coeff_of(r-i, p, q)
        return t/q(0)

def substitute_with_one_vector(gf):
    generators_all = [exp for d in gf.den_exps for exp in d]
    mu = get_non_orth_vector(generators_all)

    constant_term = 0
    for s, v, ds in gf:
        num = dot(v, mu)
        dens = [dot(d, mu) for d in ds]

        negative_exps = [d for d in dens if d < 0]
        if len(negative_exps) % 2 == 1: s = -s
        c = sum(negative_exps)

        dens = [abs(d)-1 for d in dens]
        r = len(dens)
        if r % 2 == 1: s = -s
        num -= c

        # Now we have an expression of the form
        # (s+1)^num /                => P(s)
        # s^r * Π_r((1+s)^di - 1)/s  => s^r * Q(s)
            # 1 ≤ i ≤ r; di = dens[i]+1

        p_of_s = get_coeffs_numerator(num)
        q_of_s = get_coeffs_denominator(dens)

        constant_term += s * get_coeff_of(r, p_of_s, q_of_s)

    return constant_term

def get_non_orth_vector(vectors): # none of the vectors should be 0
    dim = len(vectors[0])
    new_vector = [1]
    for d in range(2, dim+1):
        lower_dim_dot_products = [dot(new_vector, v[:d-1]) for v in vectors]
        disallowed_values = [-d/v[-1] for d, v in zip(lower_dim_dot_products, vectors) if v[-1] != 0]
        new_coordinate = ceil(max(disallowed_values) + 1).astype(int).item()
        new_vector.append(new_coordinate)
    return new_vector

def count_integer_points(poly):
    gf = polytope_gen_func(poly)
    ct = substitute_with_one_vector(gf)
    return ct

# Unit cube in first octant
mat = cdd.Matrix([[1,-1,0,0],[0,1,0,0],[1,0,-1,0],[0,0,1,0],[1,0,0,-1],[0,0,0,1]])
mat.rep_type = cdd.RepType.INEQUALITY
cube = cdd.Polyhedron(mat)
print("Cube integer points:", count_integer_points(cube))

# Tetrahedron in first octant
mat = cdd.Matrix([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, -1, -1, -1]])
mat.rep_type = cdd.RepType.INEQUALITY
tetra = cdd.Polyhedron(mat)
print("Tetrahedron integer points:", count_integer_points(tetra))