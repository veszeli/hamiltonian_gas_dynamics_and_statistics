from abc import ABC, abstractmethod
import numpy as np

class InteractionPotential(ABC):
	@abstractmethod
	def energy(self, coordinates):
		pass

	@abstractmethod
	def d_energy(self, coordinates):
		pass


class InteractionPairPotential(InteractionPotential):
	@abstractmethod
	def pair_potential(self, r1, r2):
		pass

	@abstractmethod
	def d_pair_potential(self, r1, r2):
		pass

	def energy(self, coordinates):
		n = len(coordinates)
		u = 0
		for i in range(n):
			for j in range(i, n):
				u += self.pair_potential(coordinates[i], coordinates[j])
		return u

	def d_energy(self, coordinates):
		"Partial derivative of interaction energy w.r.t. coordinates."
		n = len(coordinates)
		du = np.zeros_like(coordinates)
		for i in range(n):
			for j in np.concatenate([np.arange(0,i), np.arange(i+1,n)]):
				du[i] += self.d_pair_potential(coordinates[i],
												coordinates[j])
		return du

	


class M1Interaction(InteractionPairPotential):
	def __init__(self, r0, u0):
		"""
		$1/d^2$ like pair potential
		"""
		self.r0 = r0
		self.u0 = u0
		self.r0_sq = self.r0*self.r0
		self.b = 2*self.u0*self.r0_sq


	def pair_potential(self, r1, r2):
		d = np.linalg.norm(r1-r2)
		if d < self.r0:
			return self.u0
		else:
			return self.u0*self.r0_sq/(d*d)


	def d_pair_potential(self, r1, r2):
		d = np.linalg.norm(r1-r2)
		if d < self.r0:
			return np.zeros_like(r1)
		else:
			return -self.b*(r1-r2)/d**4


class M2Interaction(InteractionPairPotential):
	def __init__(self, u_min, u_max, d1, d2):
        """

        """
		self.u_min = u_min
		self.u_max = u_max
		self.delta_u = self.u_max - self.u_min
		self.d1 = d1
		self.d2 = d2
		self.down_slope = (self.u_min - self.u_max)/d1
		self.up_slope = -self.u_min/(self.d2-self.d1)

	def pair_potential(self, r1, r2):
		d = np.linalg.norm(r1-r2)
		if d > self.d2:
			return 0
		elif d < self.d1:
			return self.u_max + self.down_slope*d
		else:
			return self.u_min + self.up_slope*(d-self.d1)

	def d_pair_potential(self, r1, r2):
		d = np.linalg.norm(r1-r2)
		if d > self.d2:
			return 0
		elif d < self.d1:
			return self.down_slope*(r1-r2)/d
		else:
			return self.up_slope*(r1-r2)/d
