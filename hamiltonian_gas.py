import numpy as np

class HamiltonianGasSystem:
	def __init__(self, m, a, r0, u0):
		self.m = m
		self.a = a
		self.r0 = r0
		self.u0 = u0
		self.r0_sq = self.r0*self.r0
		self.b = 2*self.u0*self.r0_sq

	def energy(self, coordinates, momentums):
		return	self.kinetic_energy(momentums) +\
				self.potential_energy(coordinates) +\
				self.interaction_energy(coordinates)
	
	def kinetic_energy(self, momentums):
		return 0.5*np.sum(np.square(momentums))/self.m

	def potential_energy(self, coordinates):
		return 0.5*self.a*np.sum(np.square(coordinates))

	def pair_potential(self, r1, r2):
		d = np.linalg.norm(r1-r2)
		if d < self.r0:
			return self.u0
		else:
			return self.u0*self.r0_sq/(d*d)

	def interaction_energy(self, coordinates):
		n = len(coordinates)
		u = 0
		for i in range(n):
			for j in range(i, n):
				u += self.pair_potential(coordinates[i], coordinates[j])
		return u

	def d_kinetic(self, momentums):
		"Partial derivative of kinetic energy w.r.t. momentums."
		return momentums/self.m

	def d_potential(self, coordinates):
		"Partial derivative of potential energy w.r.t. coordinates."
		return self.a*coordinates

	def d_pair_potential(self, r1, r2):
		d = np.linalg.norm(r1-r2)
		if d < self.r0:
			return np.zeros_like(r1)
		else:
			return -self.b*(r1-r2)/d**4

	def d_interaction(self, coordinates):
		"Partial derivative of interaction energy w.r.t. coordinates."
		n = len(coordinates)
		du = np.zeros_like(coordinates)
		for i in range(n):
			for j in np.concatenate([np.arange(0,i), np.arange(i+1,n)]):
				du[i] += self.d_pair_potential(coordinates[i],
												coordinates[j])
		return du
