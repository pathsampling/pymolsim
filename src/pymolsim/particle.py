# Pure python MD code
# Based on Jutta's code
# Could be pretty slow :)
# Use only for test systems
import numpy as np

class Particle:
	"""
	Main particle class which contains info
	about a particle
	"""
	def __init__(self):
		
		self.r = np.zeros(3)
		self.v = np.zeros(3)
		self.f = np.zeros(3)
		self.type = 1
		self.mass = 1
		self.energy = 0

	def assign_pos(self, atoms):
		"""
		Takes pyscal atom objects
		"""
		pos = [atom.pos for atom in atoms]
		
		#now we have to set it in a 3XNatoms list
		self.x = np.array([p[0] for p in pos])
		self.y = np.array([p[1] for p in pos])
		self.z = np.array([p[2] for p in pos])

	def vectorize_dist(self):
		"""
		Convert the list of a distance matrix
		"""
		self.xd = np.meshgrid(self.x, self.x)[1] - np.meshgrid(self.x, self.x)[0]
		self.yd = np.meshgrid(self.y, self.y)[1] - np.meshgrid(self.y, self.y)[0]
		self.zd = np.meshgrid(self.z, self.z)[1] - np.meshgrid(self.z, self.z)[0]

		self.dist = (self.xd**2 + self.yd**2 + self.zd**2)**0.5
		


