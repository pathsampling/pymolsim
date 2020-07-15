# Pure python MD code
# Based on Jutta's code
# Could be pretty slow :)
# Use only for test systems

import numpy as np

class LJ:
	def __init__(self, epsilon, sigma):
		self.virial = 0
		self.hypervirial = 0
		self.press_kin = None
		self.energy = 0
		self.stress = np.zeros((3, 3))

		self.epsilon = epsilon
		self.sigma = sigma
		self.c1 = 0.016132*epsilon
		self.c2 = 3136.6*epsilon
		self.c3 = -68.069*epsilon
		self.c4 = -0.083312*epsilon
		self.c5 = 0.74689*epsilon

		self.rmin2 = (2.3*sigma)**2
		self.rmax2 = (2.5*sigma)**2

	def forces(self, dr):
		"""
		Calculate forces
		"""
		sigma6 = self.sigma**6
		self.virial = 0.0;
		self.hypervirial = 0.0;
		self.stress = np.zeros((3, 3))

		fi = np.zeros(3)
		fj = np.zeros(3)

		r2 = np.sum(np.array(dr)**2)

		if (r2 <= self.rmin2):
			
			r6i = 1.0/(r2*r2*r2)
			f = 24.0*self.epsilon*sigma6*r6i/r2* (2.0*sigma6*r6i - 1.0)
			g = 48.0*self.epsilon*sigma6*r6i * (13.0*sigma6*r6i - 3.5)

		elif (r2 < self.rmax2):

			r6i = 1.0/(r2*r2*r2);
			f = 2.0/r2 * (6.0*self.c2 * (sigma6*r6i)**2 + 3.0*self.c3*sigma6*r6i) - 2.0*self.c4/(self.sigma*self.sigma)
			g = 6.0*sigma6*r6i*(self.c2*26.0*sigma6*r6i + 7.0*self.c3) + 2.0*self.c4*r2/(self.sigma*self.sigma)
					
		else:
			
			f = 0.0;
			g = 0.0;
			
		self.virial += f*r2
		self.hypervirial += (g-f*r2)

		for m in range(3):
			for n in range(3):
				stress[m, n] += f*dr[m]*dr[n]

		for k in range(3):
			fi[k] -= f*dr[k]
			fj[k] += f*dr[k]

		return fi, fj


	def potential_energy(self, dr):
		"""
		Calculate pe
		"""

		energy = 0.0
		sigma6 = sigma**6

		r2 = np.sum(np.array(dr)**2)
				
		if(r2 <= self.rmin2):
			r6i = 1.0/(r2*r2*r2)
			energy += (4.0*self.epsilon*sigma6*r6i)*(sigma6*r6i - 1.0) + self.c1
		
		elif(r2 < self.rmax2):
			r6i = 1.0/(r2*r2*r2)
			energy += self.c2*(sigma6*r6i)**2 + self.c3*sigma6*r6i + self.c4*r2/(self.sigma*self.sigma) + self.c5
		
		else
			energy += 0.0

		return energy
		






