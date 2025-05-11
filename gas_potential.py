import numpy as np
from abc import ABC, abstractmethod

class OuterPotential(ABC):
	@abstractmethod
	def energy(coordinates):
		"""
		Parameters
		----------
		coordinates: numpy array, N rows, d columns

		Returns
		-------
		out: float
		"""
		pass

	def d_energy(coordinates):
		"""
		Partial derivative of potential energy w.r.t. coordinates.

		Parameters:
		-----------
		coordinates: numpy array, N rows, d columns

		Returns:
		--------
		out: numpy array, N rows, d columns
		"""
		pass


class SimpleQuadraticPotential(OuterPotential):
	def __init__(self, a):
		r"""
		Potential in the form:
		$U = a/2\sum_{i=1}^N \sum_{k=1}^d r_{ik} r_{ik}$
		Where N is the number of particles, and d is the dimension.

		Parameters
		----------
		a: float
		"""
		self.a = a

	def energy(self, coordinates):
		return self.a/2*np.sum(np.square(coordinates))

	def d_energy(self, coordinates):
		return self.a*coordinates

