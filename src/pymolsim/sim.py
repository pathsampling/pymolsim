"""
Main sim class with contains md methods
"""
import numpy as np
import os

class Sim:
	def __init__(self, atoms, boxdims):


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
		self.px = 0
		self.py = 0
		self.pz = 0
		
		self.integrator = 0
		self.potential = None
		self.thermostat = None
		self.barostat = None

		#atoms are pyscal atom objects
		pos = [atom.pos for atom in atoms]
		self.x = np.array([p[0] for p in pos])
		self.y = np.array([p[1] for p in pos])
		self.z = np.array([p[2] for p in pos])
		self.xd = None
		self.yd = None
		self.zd = None
		self.vx = np.zeros(len(atoms))
		self.vy = np.zeros(len(atoms))
		self.vz = np.zeros(len(atoms))
		self.f = np.zeros((len(atoms), len(atoms)))
		self.mass = np.ones(len(atoms))
		self.type = np.ones(len(atoms))

	def image_distance(self):

		self.xd = np.where(self.xd > self.box_2[0], self.xd-self.box[0], self.xd)
		self.xd = np.where(self.xd < -self.box_2[0], self.xd+self.box[0], self.xd)
		self.yd = np.where(self.yd > self.box_2[1], self.yd-self.box[1], self.yd)
		self.yd = np.where(self.yd < -self.box_2[1], self.yd+self.box[1], self.yd)
		self.zd = np.where(self.zd > self.box_2[2], self.zd-self.box[2], self.zd)
		self.zd = np.where(self.zd < -self.box_2[2], self.zd+self.box[2], self.zd)

	def remap(self):
		self.x = np.where(self.x > self.box[0], self.x-self.box[0], self.x)
		self.x = np.where(self.x < 0.0, self.x+self.box[0], self.x)
		self.y = np.where(self.y > self.box[1], self.y-self.box[1], self.y)
		self.y = np.where(self.y < 0.0, self.y+self.box[1], self.y)
		self.z = np.where(self.z > self.box[2], self.z-self.box[2], self.z)
		self.z = np.where(self.z < 0.0, self.z+self.box[2], self.z)

	def vectorize_dist(self):		
		"""
		The vectorise method, also need to check pbc
		"""
		self.xd = np.meshgrid(self.x, self.x)[1] - np.meshgrid(self.x, self.x)[0]
		self.yd = np.meshgrid(self.y, self.y)[1] - np.meshgrid(self.y, self.y)[0]
		self.zd = np.meshgrid(self.z, self.z)[1] - np.meshgrid(self.z, self.z)[0]
		self.image_distance()		
	
	def random_velocity(self, counts):
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
		fx, fy, fz = self.potential.forces(self.xd, self.yd, self.zd)

		#finally we should reduce the forces
		self.fx = np.sum(fx, axis=-1)
		self.fy = np.sum(fy, axis=-1)
		self.fz = np.sum(fz, axis=-1)

	def potential_energy(self):

		energy = self.potential.potential_energy(self.xd, self.yd, self.zd)
		pe = np.sum(energy)/2
		#we need to divide by 2 to remove the overcounting
		return pe


	def kinetic_energy(self):
		ke = np.sum(self.mass*(self.vx**2 + self.vy**2 + self.vz**2))
		return ke/2

	def total_energy(self):
		pe = self.potential_energy()
		ke = self.kinetic_energy()
		return pe+ke

	def total_momentum(self):
		self.px = np.sum(self.mass*self.vx)
		self.py = np.sum(self.mass*self.vy)
		self.pz = np.sum(self.mass*self.vz)

	def rescale_velocities(self):
		ke = self.kinetic_energy()
		temp = 2.0*ke/(3*self.nparticles - 3.0)
		temp_aim = 1.0/self.beta
		c = np.sqrt(temp_aim/temp)

		self.vx = self.vx*c
		self.vy = self.vy*c
		self.vz = self.vz*c

	@property
	def temperature(self):
		temp = 2.0*self.kinetic_energy()/(3*self.nparticles - 3.0)
		return temp


	def start(self):
		
		self.vx = self.random_velocity(counts=self.nparticles)
		self.vy = self.random_velocity(counts=self.nparticles)
		self.vz = self.random_velocity(counts=self.nparticles)

		#find momenta
		self.total_momentum()
		ke = self.kinetic_energy()

		self.vx -= self.px/self.nparticles/self.mass
		self.vy -= self.py/self.nparticles/self.mass
		self.vz -= self.pz/self.nparticles/self.mass

		#assign thermostats
		if self.thermostat is None:
			self.run = self.md_verlet
			print("running MD Verlet")
		elif self.thermostat.name == "rescale":
			self.run = self.md_rescale
			print("running with temp rescaling")
		elif self.thermostat.name == "andersen":
			self.run = self.md_andersen
			print("running with Andersen thermostat")
		elif self.thermostat.name == "langevin":
			self.run = self.md_langevin
			print("running with Langevin thermostat")
		
		#pe = self.potential_energy()
	def md_verlet(self):
		self.propagate_momenta_half()
		self.propagate_position_half()
		self.propagate_position_half()
		self.forces()
		self.propagate_momenta_half()
		self.remap()

	def md_rescale(self):
		self.md_verlet()
		self.rescale_velocities()

	def md_langevin(self):
		self.langevin_thermo()
		self.md_verlet()
		self.langevin_thermo()
		#self.remap()

	def md_andersen(self):
		self.md_verlet()
		rands = np.random.rand(self.nparticles)
		
		pvx = np.sqrt(1.0/(self.mass*self.beta))*np.random.normal()
		pvy = np.sqrt(1.0/(self.mass*self.beta))*np.random.normal()
		pvz = np.sqrt(1.0/(self.mass*self.beta))*np.random.normal()
		
		self.vx = np.where(rands < self.thermostat.anu*self.dt, pvx, self.vx)
		self.vy = np.where(rands < self.thermostat.anu*self.dt, pvy, self.vy)
		self.vz = np.where(rands < self.thermostat.anu*self.dt, pvz, self.vz)
		#self.remap()

	def langevin_thermo(self):
		self.vx = self.thermostat.lc1*self.vx + self.thermostat.lc2/np.sqrt(self.mass)*np.random.normal(size=self.nparticles)
		self.vy = self.thermostat.lc1*self.vy + self.thermostat.lc2/np.sqrt(self.mass)*np.random.normal(size=self.nparticles)
		self.vz = self.thermostat.lc1*self.vz + self.thermostat.lc2/np.sqrt(self.mass)*np.random.normal(size=self.nparticles)

	def propagate_momenta_half(self):
		self.vx += 0.5*self.dt*self.fx/self.mass
		self.vy += 0.5*self.dt*self.fy/self.mass
		self.vz += 0.5*self.dt*self.fz/self.mass

	def propagate_position_half(self):
		self.x += 0.5*self.dt*self.vx
		self.y += 0.5*self.dt*self.vy
		self.z += 0.5*self.dt*self.vz		
	
	
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
			fout.write("%f %f %f "%(self.x[i], self.y[i], self.z[i]))
			fout.write("%f %f %f"%(self.vx[i], self.vy[i], self.vz[i]))
			fout.write("\n")
		
		fout.close()			


