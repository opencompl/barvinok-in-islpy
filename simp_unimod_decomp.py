import olll
import numpy as np

# For simplicial cones
class Cone():
    def __init__(self, rays):
        self.rays = rays # [d, d] matrix for simplicial cones
    
    def get_index(self):
        """Volume of parallelepiped"""
        det = np.linalg.det(self.rays)
        return (0.5 * det)

def unimodular_decomp(cone):
    ind = cone.get_index()
    if ind == 1: return [cone]
    else:
        cones = []
        for i in range(len(cone.rays)):
            w = "integer point lying in {∑ai ui | ai ≤ ind}"
            ki = Cone("u[0], ..., u[i-1], w, u[i+1], ..., u[-1]")
            cones.append(ki)
        final_nested = map(unimodular_decomp, cones)
        final = [cone for decomp in final_nested for cone in decomp]
        return final