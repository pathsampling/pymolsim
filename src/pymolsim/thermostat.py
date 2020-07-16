"""
Thermostats for use in 
pymolsim module
"""
import numpy as np

class Rescale:
	def __init__(self):
		self.name = "rescale"

class Andersen:
	def __init__(self, frequency=0.0):
		self.name = "andersen"
		self.anu = frequency

class Langevin:
	def __init__(self, gamma=0.0):
		self.name = "langevin"
		self.lgamma = gamma

	def start(self, dt, beta):
		self.lc1 = np.exp(-self.lgamma*dt/2.0)
		self.lc2 = np.sqrt((1.0-(self.lc1**2))/beta)




