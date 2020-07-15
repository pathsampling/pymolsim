# Pure python MD code
# Based on Jutta's code
# Could be pretty slow :)
# Use only for test systems
import numpy as np

class Barostat:
	def __init__(self,):
		self.isotropic = True
		self.pv = None
		self.pmass = None
		self.lgamma = None

	def start(self, dt, beta):
		self.lc1 = np.exp(-self.lgamma*dt/2.0)
		lc2 = np.sqrt((1.0-(self.lc1**2))/beta)
		self.lc2 = lc2/sqrt(self.pmass)


