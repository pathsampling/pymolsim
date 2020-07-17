"""
Main sim class with contains md methods
"""
import numpy as np
import os

class Sim:
	def __init__(self, atoms, boxdims):

		#note that you should still give a box in all three dimensions to trick
		self.nparticles = len(atoms)
		self.boxdims = boxdims
		self.box = np.array([boxdims[0][1]-boxdims[0][0], boxdims[1][1]-boxdims[1][0], boxdims[2][1]-boxdims[2][0]])
		self.box_2 = np.array([(boxdims[0][1]-boxdims[0][0])/2, (boxdims[1][1]-boxdims[1][0])/2, (boxdims[2][1]-boxdims[2][0])/2])
		self.box_relative = np.array([1.0, self.box[1]/self.box[0], self.box[2]/self.box[0]])
		self.volume = self.box[0]*self.box[1]*self.box[2]
		self.rho = self.nparticles/self.volume
		self.trajfile = os.path.join(os.getcwd(), "traj.dat")

		self.beta = None
		self.pressure = None
		self.dt = 0.0025
		self.time = 0
		self.tot_energy = 0
		self.pe = 0
		self.ke = 0

		self.p = np.zeros(3)
		#self.px = 0
		#self.py = 0
		#self.pz = 0
		
		self.integrator = 0
		self.potential = None
		self.thermostat = None
		self.barostat = None

		#atoms are pyscal atom objects
		pos = [atom.pos for atom in atoms]
		self.x = np.array([ np.array([p[0] for p in pos]),
				            np.array([p[1] for p in pos]),
				            np.array([p[2] for p in pos]) ])

		#self.x = np.array([p[0] for p in pos])
		#self.y = np.array([p[1] for p in pos])
		#self.z = np.array([p[2] for p in pos])
		self.xd = None
		#self.xd = None
		#self.yd = None
		#self.zd = None
		self.v = None
		
		#self.vx = np.zeros(len(atoms))
		#self.vy = np.zeros(len(atoms))
		#self.vz = np.zeros(len(atoms))
		self.f = np.zeros((3, len(atoms), len(atoms)))
		self.mass = np.ones(len(atoms))
		self.type = np.ones(len(atoms))
		self.dim = 3

	def image_distance(self):
		for i in range(self.dim):
			self.xd[i] = np.where(self.xd[i] > self.box_2[i], self.xd[i]-self.box[i], self.xd[i])
			self.xd[i] = np.where(self.xd[i] < -self.box_2[i], self.xd[i]+self.box[i], self.xd[i])

	def remap(self):
		for i in range(self.dim):
			self.x[i] = np.where(self.x[i] > self.box[i], self.x[i]-self.box[i], self.x[i])
			self.x[i] = np.where(self.x[i] < 0.0, self.x[i]+self.box[i], self.x[i])

	def vectorize_dist(self):		
		"""
		The vectorise method, also need to check pbc
		"""
		xdum = []
		for i in range(self.dim):
			xd = np.meshgrid(self.x[i], self.x[i])[1] - np.meshgrid(self.x[i], self.x[i])[0]
			xdum.append(xd)
		self.xd = np.array(xdum)
					
		self.image_distance()		
	
	def random_velocity(self, counts=1):
		"""
		Generate random velocities
		"""
		temperature = 1.0/self.beta
		vel = np.sqrt(temperature/self.mass)*np.random.normal(size=counts)
		return vel

	def forces(self):
		"""
		Clauclate forces of the system
		"""
		#first step to vectorize dist
		self.vectorize_dist()

		#now calculate forces
		fx = self.potential.forces(self.xd, dim=self.dim)

		#finally we should reduce the forces
		fdum = []
		for i in range(self.dim):
			fdum.append(np.sum(fx[i], axis=-1))
		self.f = np.array(fdum)

	def potential_energy(self):

		energy = self.potential.potential_energy(self.xd, dim = self.dim)
		pe = np.sum(energy)/2
		#we need to divide by 2 to remove the doublecounting
		return pe


	def kinetic_energy(self):
		ke = np.sum(self.mass*(self.v[0]**2 + self.v[1]**2 + self.v[2]**2))
		return ke

	def total_energy(self):
		pe = self.potential_energy()
		ke = self.kinetic_energy()
		return pe+ke

	def total_momentum(self):
		for i in range(self.dim):
			self.p[i] = np.sum(self.mass*self.v[i])

	def rescale_velocities(self):
		ke = self.kinetic_energy()
		temp = 2.0*ke/(self.dim*self.nparticles - self.dim)
		temp_aim = 1.0/self.beta
		c = np.sqrt(temp_aim/temp)

		for i in range(self.dim):
			self.v[i] = self.v[i]*c

	@property
	def temperature(self):
		temp = 2.0*self.kinetic_energy()/(self.dim*self.nparticles - self.dim)
		return temp


	def start(self):
		
		dumv = []
		for i in range(self.dim):
			dumv.append(self.random_velocity(counts=self.nparticles))
		self.v = np.array(dumv)

		#find momenta
		self.total_momentum()
		ke = self.kinetic_energy()

		for i in range(self.dim):
			self.v[i] -= self.p[i]/self.nparticles/self.mass

		#assign thermostats
		if self.thermostat is None:
			self.run = self.md_verlet
		elif self.thermostat.name == "rescale":
			self.run = self.md_rescale
		elif self.thermostat.name == "andersen":
			self.run = self.md_andersen
		elif self.thermostat.name == "langevin":
			self.run = self.md_langevin
		
		#pe = self.potential_energy()
	def md_verlet(self):
		self.propagate_momenta_half()
		self.propagate_position_half()
		self.propagate_position_half()
		self.forces()
		self.propagate_momenta_half()
		self.remap()

	def md_rescale(self):
		self.propagate_momenta_half()
		self.propagate_position_half()
		self.propagate_position_half()
		self.forces()
		self.propagate_momenta_half()
		self.remap()
		self.rescale_velocities()

	def md_langevin(self):
		self.langevin_thermo()
		self.md_verlet()
		self.langevin_thermo()
		self.remap()

	def md_andersen(self):
		self.md_verlet()
		rands = np.random.rand(self.nparticles)
		for i in range(self.dim):
			self.v[i] = np.where(rands < self.thermostat.anu*self.dt, np.sqrt(1.0/(self.mass*self.beta))*np.random.normal(), self.v[i])
		self.remap()

	def langevin_thermo(self):
		for i in range(self.dim):
			self.v[i] = self.thermostat.lc1*self.v[i] + self.thermostat.lc2/np.sqrt(self.mass)*np.random.normal(size=self.nparticles)

	def propagate_momenta_half(self):
		for i in range(self.dim):
			self.v[i] += 0.5*self.dt*self.f[i]/self.mass

	def propagate_position_half(self):
		for i in range(self.dim):
			self.x[i] += 0.5*self.dt*self.v[i]
	
	
	def dump(self, step):
		"""
		Cant vectorize this :)
		"""
		
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

		for i, particle in enumerate(self.x):
			fout.write("%d %d %d "%(i+1, self.type[i], self.mass[i]))
			fout.write("%f %f %f "%(self.x[0][i], self.x[1][i], self.x[2][i]))
			fout.write("%f %f %f"%(self.v[0][i], self.v[1][i], self.v[2][i]))
			fout.write("\n")
		
		fout.close()			


