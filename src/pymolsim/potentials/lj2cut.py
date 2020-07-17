import numpy as np

class LJ:
	def __init__(self, epsilon, sigma):

		self.epsilon = epsilon
		self.sigma = sigma
		self.c1 = 0.016132*epsilon
		self.c2 = 3136.6*epsilon
		self.c3 = -68.069*epsilon
		self.c4 = -0.083312*epsilon
		self.c5 = 0.74689*epsilon

		self.rmin2 = (2.3*sigma)**2
		self.rmax2 = (2.5*sigma)**2

	def forces(self, xd, dim=3):
		"""
		Vectorized force call
		input is an natomsxnatoms  arrays
		"""
		#We need to vectorize the function call - so we need to add some subfunctions
		r2 = xd[0]**2
		for i in range(1, dim):
			r2 += (xd[i]**2)
		r6i = np.where(np.abs(r2)>0, 1/(r2*r2*r2), r2)

		def _force_cut1(r2):
			sigma6 = self.sigma**6
			r6j = np.where(np.abs(r2)>0, r6i/r2, r6i)
			f = 24.0*self.epsilon*sigma6*r6j*(2.0*sigma6*r6i - 1.0)
			return f

		def _force_cut2(r2):
			sigma6 = self.sigma**6
			r2i = np.where(np.abs(r2)>0, 2.0/r2, r2)			
			f = r2i * (6.0*self.c2 * (sigma6*r6i)**2 + 3.0*self.c3*sigma6*r6i) - 2.0*self.c4/(self.sigma*self.sigma)
			return f

		fnull = np.zeros_like(r2)
		fcut1 = np.where(r2 <= self.rmin2, _force_cut1(r2), fnull)
		fcut2 = np.where((r2 > self.rmin2) & (r2 < self.rmax2), _force_cut2(r2), fnull)
		f = fcut1 + fcut2

		#convert forces to its components
		fx = []
		for i in range(dim):
			fx.append(f*xd[i])
		
		return np.array(fx)
		

	def potential_energy(self, xd, dim=3):
		"""
		Calculate pe : vectorized
		"""
		r2 = xd[0]**2
		for i in range(1, dim):
			r2 += (xd[i]**2)

		def _pe_cut1(r2):
			sigma6 = sigma**6
			r6i = 1.0/(r2*r2*r2)
			energy = (4.0*self.epsilon*sigma6*r6i)*(sigma6*r6i - 1.0) + self.c1
			return energy

		def _pe_cut2(r2):
			sigma6 = sigma**6				
			r6i = 1.0/(r2*r2*r2)
			energy = self.c2*(sigma6*r6i)**2 + self.c3*sigma6*r6i + self.c4*r2/(self.sigma*self.sigma) + self.c5
			return energy
		
		energy = np.zeros_like(r2)
		energy = np.where(r2 <= self.rmin2, _pe_cut1(r2), energy)
		energy = np.where((r2 > self.rmin2) & (r2 < self.rmax2), _pe_cut2(r2), energy)

		return energy