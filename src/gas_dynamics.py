import numpy as np

class GasDynamics:
    def __init__(self, d_hamiltonian_dr, d_hamiltonian_dp, half_box_size):
        self.f = d_hamiltonian_dp
        self.g = lambda x: -d_hamiltonian_dr(x)
        self.d = half_box_size

    def sie(self, r, p, dt):
        "Semi implicit Euler method"
        p2 = p + self.g(r)*dt
        r2 = r + self.f(p2)*dt
        return r2, p2

    def walls(self, r,p):
        "Effect of walls: particles bounces back"
        num_particles = len(r)
        dimension = len(r[0])
        for i in range(num_particles):
            for k in range(dimension):
                if r[i][k] > self.d:
                    r[i][k] = self.d
                    p[i][k] *= -1
                if r[i][k] < -self.d:
                    r[i][k] = -self.d
                    p[i][k] *= -1

                    

