"""
This would implement a simple double well potential
as a toy model for TIS sampling.

Going through a revamp to constitute a general potential

15Feb20 - Major restructuring to make the module more generic
functions. For example, MD in 2, 3 or 1 dimensions.
Cancelled - we stick to two dimensions

But -> the output files and so on would depend on the potential
themselves - we might not need all the information.

We define singleparticle and multiparticle classes. This in turn
would contain different potentials and also options to write and
read trajectories.

"""
import sys
import os
import numpy as np
from numba import njit


class Atom:
    """
    Basic data structure to store atom details
    Mass is 1 for now
    """
    def __init__(self, x, y, vx, vy):
        """
        r, v, f and m position, velocity, force and mass
        """
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
        self.fx = None
        self.fy = None
        self.mass = 1
        self.type = 1
        self.neighbors = None
        self.diffxs = None
        self.diffys = None


class MdRoutine(object):

    def __init__(self,  potential, atoms=None, gamma=2.5, beta=2, timestep=0.05, box=[-6., 6., -4., 4.], overlay_potential=None):
        """
        Potential is to be converted to a list?
        """
        self.atoms = atoms
        self.potential = potential
        self.dt = timestep
        self.boxx1 = box[0]
        self.boxx2 = box[1]
        self.boxy1 = box[2]
        self.boxy2 = box[3]
        self.gamma = gamma
        self.beta = beta
        self.lc1 = np.exp(-self.gamma*self.dt/2)
        self.lc2 = np.sqrt((1 - self.lc1**2)/self.beta)
        self.delangevin = 0.0
        #if more than one potential is to be used, overlay_potential can be set
        self.overlay_potential = overlay_potential

    def update_forces(self):
        """
        Update forces on atoms
        """
        for atom in self.atoms:
            energy, fx, fy = self.potential(atom.x, atom.y)
            atom.fx = fx
            atom.fy = fy


    def propagate_velocity_half(self):
        """
        propagate velocity by half time step
        """
        for atom in self.atoms:
            atom.vx = atom.vx + 0.5*(atom.fx/atom.mass)*self.dt
            atom.vy = atom.vy + 0.5*(atom.fy/atom.mass)*self.dt

    def propagate_position_half(self):
        """
        Propagate position by half time step
        """
        for atom in self.atoms:
            atom.x = atom.x + 0.5*atom.vx*self.dt
            atom.y = atom.y + 0.5*atom.vy*self.dt

    def md_verlet(self):
        """
        Do a verlet half step
        """
        self.propagate_velocity_half()
        self.propagate_position_half()
        self.propagate_position_half()
        #find forces
        self.update_forces()
        self.propagate_velocity_half()

    def md_langevin(self):
        """
        Do a langevin md step
        """
        self.langevin_thermo()
        self.md_verlet()
        self.langevin_thermo()


    def langevin_thermo(self):
        """
        Langevin half step integrator-
        refer Jutta's code for advanced atomistic exercise
        Also - Bussi, Parrinello, PRE 75, 056707 (2007)
        """
        # update velocities and account for langevin contribution
        # to the kinetic energy (to be substracted later)
        for atom in self.atoms:
            self.delangevin += 0.5*atom.mass*atom.vx**2
            atom.vx = self.lc1*atom.vx + self.lc2*np.random.normal()
            self.delangevin -= 0.5*atom.mass*atom.vx**2

            self.delangevin += 0.5*atom.mass*atom.vy**2
            atom.vy = self.lc1*atom.vy + self.lc2*np.random.normal()
            self.delangevin -= 0.5*atom.mass*atom.vy**2
