# Pure python MD code
# Based on Jutta's code
# Could be pretty slow :)
# Use only for test systems

import numpy as np

class Thermostat:
	def __init__(self, ):
		self.lgamma = None
		self.lc1 = None
		self.lc2 = None
		self.anu = None
		self.nhlgamma = None
		self.nhlmu = None
		self.nhlc1 = None
		self.nhlc2 = None
		self.nhlxi = None
		self.dElangevin = 0

	def start(self, dt, beta, integrator):

		self.lc1 = np.exp(-self.lgamma*dt/2.0)
		self.lc2 = np.sqrt((1.0-(self.lc1**2))/beta)

		#self.nhlc1 = exp(-self.nhlgamma*dt/2.0);
		#self.nhlc2 = sqrt((1.0-(self.nhlc1**2))*2.0/beta/self.nhlmu);        

		#self.nhlxi = random_velocity(beta,self.nhlmu)

	def random_velocity(self, beta, mass):
		"""
		Generate random velocities
		"""
		temperature = 1.0/beta
		vel = np.sqrt(temperature/mass)*np.random.normal()
		return vel





