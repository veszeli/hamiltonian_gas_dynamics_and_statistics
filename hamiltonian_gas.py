import numpy as np


class HamiltonianGasSystem:
	def __init__(self, m, outer_potential, interaction):
		"""
		outer_potenetial: OuterPotential
		interaction: InteractionPotential
		"""
		self.m = m
		self.outer_potential = outer_potential
		self.interaction = interaction


	def energy(self, coordinates, momentums):
		return	self.kinetic_energy(momentums) +\
				self.outer_potential.energy(coordinates) +\
				self.interaction.energy(coordinates)
	
	def kinetic_energy(self, momentums):
		return 0.5*np.sum(np.square(momentums))/self.m


	def d_kinetic(self, momentums):
		"Partial derivative of kinetic energy w.r.t. momentums."
		return momentums/self.m


	def d_hamiltonian_d_coordinates(self, coordinates):
		return self.outer_potential.d_energy(coordinates) +\
				 self.interaction.d_energy(coordinates)


	def d_hamiltonian_d_momentums(self, momentums):
		return self.d_kinetic(momentums)



