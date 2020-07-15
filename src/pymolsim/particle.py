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
		