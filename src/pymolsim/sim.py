# Pure python MD code
# Based on Jutta's code
# Could be pretty slow :)
# Use only for test systems
import numpy as np
import os

class Sim:
	def __init__(self, particles, boxdims):
		
		self.nparticles = len(particles)
		self.boxdims = boxdims
		self.box = np.array([boxdims[0][1]-boxdims[0][0], boxdims[1][1]-boxdims[1][0], boxdims[2][1]-boxdims[2][0]])
		self.box_2 = np.array([(boxdims[0][1]-boxdims[0][0])/2, (boxdims[1][1]-boxdims[1][0])/2, (boxdims[2][1]-boxdims[2][0])/2])
		self.box_relative = np.array([1.0, self.box[1]/self.box[0], self.box[2]/self.box[0]])
		self.volume = self.box[0]*self.box[1]*self.box[2]
		self.rho = self.nparticles/self.volume
		self.trajfile = os.path.join(os.getcwd(), "traj.dat")
		self.particles = particles

		self.E0 = 0
		self.beta = None
		self.pressure = None
		self.dt = 0.0025
		self.time = 0
		self.tot_energy = 0
		self.pe = 0
		self.ke = 0
		self.p = np.zeros(3)

		self.Tinst = 0
		self.Pinst = 0
		self.Pinst_1 = 0
		self.virial = 0
		self.hypervirial = 0
		self.press_kin = 0
		self.stress = np.zeros((3, 3))

		self.integrator = 0
		self.potential = None
		self.thermostat = None
		self.barostat = None
		self.nthreads = 1
		
	
	def random_velocity(self, mass):
		"""
		Generate random velocities
		"""
		temperature = 1.0/self.beta
		vel = np.sqrt(temperature/mass)*np.random.normal()
		return vel

	def forces(self):
		"""
		Clauclate forces of the system
		"""
		for particle in self.particles:
			particle.f = np.zeros(3)

		#now this is the main loop - daskify later
		calcs = []
		for i in range(self.nparticles - 1):
			for j in range(i+1, self.nparticles):
				self.call_force(i, j)

		self.stress += self.potential.stress/self.volume
		self.virial = self.potential.virial/(3*self.volume)
		self.hypervirial = self.potential.hypervirial/9

				
		
	def call_force(self, i, j):
		"""
		Sub force call
		"""
		dr = self.image_distance(i, j)
		fi, fj = self.potential.forces(dr)
		#These are vector operations
		self.particles[i].f += fi
		self.particles[j].f += fj

	def potential_energy(self):
		pe = 0
		for i in range(self.nparticles - 1):
			for j in range(i+1, self.nparticles):
				pe += self.call_potential_energy(i, j)


	def call_potential_energy(self, i, j):
		dr = self.image_distance(i, j)
		return self.potential.potential_energy(dr)


	def kinetic_energy(self):
		ke = 0
		for particle in self.particles:
			ke += particle.mass*np.sum(particle.v**2)
		return ke

	def total_energy(self):
		pe = self.potential_energy()
		ke = self.kinetic_energy()
		return pe+ke

	def total_momentum(self):
		self.p = np.zeros(3)
		for particle in self.particles:
			self.p += particle.v*particle.mass

	def rescale_velocities(self):
		ke = self.kinetic_energy()
		temp = 2.0*ke/(3*self.nparticles - 3.0)
		temp_aim = 1.0/self.beta
		c = np.sqrt(temp_aim/temp)

		for particle in self.particles:
			particle.v *= c

	def image_distance(self, i, j):
		dr = self.particles[i].r - self.particles[j].r
		for i in range(3):
			if (dr[i] > self.box_2[i]):
				dr[i] -= self.box[i]
			elif (dr[i] < -self.box_2[i]):
				dr[i] += self.box[i]
		return dr

	def remap(self):
		for particle in self.particles:
			for i in range(3):
				if (particle.r[i] > self.box[i]):
					particle.r[i] -= self.box[i]
				if (particle.r[i] < 0.0):
					particle.r[i] += self.box[i]

	def start(self):
		for particle in self.particles:
			for i in range(3):
				particle.v[i] = self.random_velocity(particle.mass)
				self.p[i] += particle.v[i]*particle.mass

		ke = self.kinetic_energy()

		self.p /= self.nparticles

		for particle in self.particles:
			particle.v = self.p/particle.mass

		
		#pe = self.potential_energy()
	def md_verlet(self):
		self.propagate_momenta_half()
		self.propagate_position_half()
		self.propagate_position_half()
		self.forces()
		self.propagate_momenta_half()

	def md_langevin(self):
		self.langevin_thermo()
		self.md_verlet()
		self.langevin_thermo()

	def md_andersen(self):
		self.md_verlet()
		for particle in self.particles:
			if np.random.rand() < self.thermostat.anu*self.dt:
				for i in range(3):
					particle.v[i] = np.sqrt(1.0/(particle.mass*self.beta))*np.random.normal()
	
	def md_nosehooverlangevin_NVT(self):
		Nf = 3*self.nparticles - 3
		self.thermostat.nhlxi = self.thermostat.nhlc1*self.thermostat.nhlxi + self.thermostat.nhlc2*np.random.normal();
		self.propagate_momenta_xi()
		self.propagate_momenta_half()
		self.propagate_position_half()

		ke = self.kinetic_energy()
		self.thermostat.nhlxi = self.thermostat.nhlxi + self.dt*(2.0*ke - (Nf/self.beta))/self.thermostat.nhlmu
		self.propagate_position_half()
		self.forces()
		self.propagate_momenta_half()
		self.propagate_momenta_xi()
		self.thermostat.nhlxi = self.thermostat.nhlc1*self.thermostat.nhlxi + self.thermostat.nhlc2*np.random.normal()
	

	def langevin_thermo(self):
		for particle in self.particles:
			for i in range(3):
				self.thermostat.dElangevin += 0.5*particle.mass*particle.v[i]**2;
				particle.v[i] = self.thermostat.lc1*particle.v[i] + self.thermostat.lc2/np.sqrt(particle.mass)*np.random.normal()
				self.thermostat.dElangevin -= 0.5*particle.mass*particle.v[i]**2

	def langevin_baro(self):
		self.barostat.pv = self.barostat.lc1*self.barostat.pv + self.barostat.lc2*np.random.normal()

	def propagate_momenta_half(self):
		for particle in self.particles:
			particle.v += 0.5*self.dt*particle.f/particle.mass

	def propagate_position_half(self):
		for particle in self.particles:
			particle.r += 0.5*self.dt*particle.v		
	
	def propagate_momenta_xi(self):
		c = np.exp(-self.thermostat.nhlxi*self.dat/2)
		for particle in self.particles:
			particle.v = particle.v *c

	def langevin_xi(self):
		self.thermostat.nhlxi = self.thermostat.nhlc1*self.thermostat.nhlxi + self.thermostat.nhlc2*np.random.normal()
	
	#def propagate_position_half_scale(self):
	#def rescale_position_momenta(self):
	
	def dump(self, step):
		
		if (step == 0):
			fout = open(self.trajfile, "w")
		else:
			fout = open(self.trajfile, "a")	
		
		fout.write("ITEM: TIMESTEP\n")
		fout.write("%d\n"%step)
		fout.write("ITEM: NUMBER OF ATOMS\n")
		fout.write("%d\n"%self.nparticles)
		fout.write("ITEM: BOX BOUNDS pp pp pp\n")
		fout.write("0 %f\n"%self.box[0])
		fout.write("0 %f\n"%self.box[1])
		fout.write("0 %f\n"%self.box[2])
		fout.write("ITEM: ATOMS id type mass x y z vx vy vz\n")

		for i, particle in enumerate(self.particles):
			fout.write("%d %d %d "%(i+1, particle.type, particle.mass))
			for j in range(3):
				fout.write("%f "%particle.r[j])
			for j in range(3):
				fout.write("%f "%particle.v[j])
			fout.write("\n")
		
		fout.close()			


