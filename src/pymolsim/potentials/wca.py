"""
Sample system of N repulsive particles in 
Weeks-Chandler-Andersen potential. 

See ref J. Chem. Phys., Vol. 118, No. 17, 1 May 2003
"""

import numpy as np

class WCA:
	def __init__(self, epsilon=1.0, sigma=1.0, h=15, w=0.5):

		self.epsilon = epsilon
		self.sigma = sigma
		self.h = h
		self.w = w
		self.rcutsq = (self.sigma*(2)**(1/6))**2


	def forces(self, xd, types, dim=2):
		"""
		Vectorized force call
		input is an natomsxnatoms  arrays
		"""
		#We need to vectorize the function call - so we need to add some subfunctions
		r2 = xd[0]**2
		for i in range(1, dim):
			r2 += (xd[i]**2)
		
		sigma6 = self.sigma**6
		sigma12 = sigma6**2
		r6i = np.where(np.abs(r2)>0, 1/(r2*r2*r2), r2)
		r6j = np.where(np.abs(r2)>0, r6i/r2, r6i)

		def _force_cut(r2):
			f = 24.0*self.epsilon*sigma6*r6j*(2.0*sigma6*r6i - 1.0)
			return f

		fnull = np.zeros_like(r2)
		f = np.where(r2 <= self.rcutsq, _force_cut(r2), fnull)
		
		#convert forces to its components
		fx = []
		for i in range(dim):
			fx.append(f*xd[i])

		#now we have to add the extra force on pairs of atoms
		#find where the atoms are of type 2
		t2atoms = np.where(np.array(types)==2)[0]
		#get the distance between t2atoms
		t2dist = r2[t2atoms[0], t2atoms[1]]
		#great, now calculate the force between the two atoms
		r0 = self.cutsq**0.5
		t2f = (4*self.h/self.w**2)*(t2dist - r0 - self.w)*(1 - ((t2dist - r0 - self.w)/self.w)**2)
		#add the forces
		for i in range(dim):
			fx[i][t2atoms[0], t2atoms[1]] += t2f*xd[i]
			fx[i][t2atoms[1], t2atoms[0]] -= t2f*xd[i]
		
		return np.array(fx)
		

	def potential_energy(self, xd, types, dim=2):
		"""
		Calculate pe : vectorized
		"""
		r2 = xd[0]**2
		for i in range(1, dim):
			r2 += (xd[i]**2)

		def _pe_cut1(r2):
			sigma6 = sigma**6
			r6i = 1.0/(r2*r2*r2)
			energy = (4.0*self.epsilon*sigma6*r6i)*(sigma6*r6i - 1.0) + self.epsilon
			return energy
		
		energy = np.zeros_like(r2)
		energy = np.where(r2 <= self.rcutsq, _pe_cut1(r2), energy)

		#find where the atoms are of type 2
		t2atoms = np.where(np.array(types)==2)[0]
		#get the distance between t2atoms
		t2dist = r2[t2atoms[0], t2atoms[1]]
		#great, now calculate the force between the two atoms
		r0 = self.cutsq**0.5
		t2e =  self.h*(1- ((t2dist - r0 - self.w)/self.w)**2)**2
		#add the forces
		energy[t2atoms[0], t2atoms[1]] += t2e
		energy[t2atoms[1], t2atoms[0]] += t2e

		return energy