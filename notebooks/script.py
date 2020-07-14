import pyscal.core as pc
from pymolsim.cmolsim import Sim, Thermostat, Barostat, Average, Particle, Potential
import pyscal.crystal_structures as pcs

sys = pc.System()
sys.read_inputfile("conf.dump", customkeys=["mass", "vx", "vy", "vz"])
atoms = sys.atoms
box = sys.box
#convert to particles
particles = []
for atom in atoms:
    p = Particle()
    p.mass = 1.00;
    p.r = atom.pos;
    p.v = [float(atom.custom["vx"]), float(atom.custom["vy"]), float(atom.custom["vz"])]
    p.type = atom.type
    particles.append(p)
sim = Sim()
sim.particles = particles
sim.nparticles = len(particles)
sim.box = [sys.box[0][1]-sys.box[0][0], sys.box[1][1]-sys.box[1][0], sys.box[2][1]-sys.box[2][0]]
sim.box_2 = [(sys.box[0][1]-sys.box[0][0])/2, (sys.box[1][1]-sys.box[1][0])/2, (sys.box[2][1]-sys.box[2][0])/2]
sim.box_relative = [1.0, sim.box[1]/sim.box[0], sim.box[2]/sim.box[0]]
sim.volume = sim.box[0]*sim.box[1]*sim.box[2]
sim.rho = sim.nparticles/sim.volume

sim.integrator = 0
sim.beta = 0.5
sim.pressure = 0
sim.dt = 0.0025

pot = Potential()
pot.epsilon = 1.0
pot.sigma = 1.0
pot.init()
sim.dt = 0.0025
sim.potential = pot
sim.init()
sim.dt = 0.0025
sim.forces()
sim.md_step()
print(sim.dt)
