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
from pathsampling.potentialsmd import Atom

class MultiParticle:
    """
    Class contains multi particle potentials and methods to write/read trajectories
    """
    def __init__(self):
        pass

    def wca(self, r, **kwargs):
        """
        Week-Chandler-Andersen potential. See JCP 54, 5237, 1971
        sigma, epsilon, rcut are the parameters

        r - is a single component distance

        2D potential
        """
        sigma = kwargs.get('sigma', 1)
        epsilon = kwargs.get('epsilon', 1)
        rcut = kwargs.get('rcut', 1.1225)

        potential = 0
        force = 0

        if (r <= rcut):
            #only if within cutoff distance
            t1 = (r/sigma)**-6
            t2 = t1**2
            potential = 4*epsilon*(t2 - t1) + epsilon

            t1_dx = -6*t1*(r/sigma)**-1
            t2_dx = 2*t1*t1_dx
            force = -(4*epsilon*(t2_dx - t1_dx))

        return potential, force

    def double_well(self, r, **kwargs):
        """
        Double well potential
        """
        h = kwargs.get('h', 15)
        w = kwargs.get('w', 0.5)
        rcut = kwargs.get('rcut', 1.1225)

        t1 = (1 - ((r - rcut - w)**2/w**2))
        potential = h*t1**2

        t1_dx = h*2*t1*-1*(1/w**2)*2*(r-rcut-w)
        force = -t1_dx

        return potential, force

    def get_distance(self, diffx, diffy, boxx, boxy):
        """
        Check periodic boundary conditions
        """
        if (diffx > boxx/2):
            diffx -= boxx
        if (diffx < -boxx/2):
            diffx += boxx
        if (diffy > boxy/2):
            diffy -= boxy
        if (diffy < -boxy/2):
            diffy += boxy

        #claculate distance
        dist = np.sqrt(diffx*diffx + diffy*diffy)
        return dist, diffx, diffy



    def find_neighbors(self, rcut, atoms, box):
        """
        Fast method to calculate interatomic distance
        """
        boxx = box[0][1] - box[0][0]
        boxy = box[1][1] - box[1][0]

        for i in range(len(atoms)):
            atoms[i].neighbors = []
            atoms[i].diffxs = []
            atoms[i].diffys = []

        for i in range(len(atoms)):
            for j in range(i+1, len(atoms)):
                diffx = atoms[i].x - atoms[j].x
                diffy = atoms[i].y - atoms[j].y
                dist, diffx, diffy = self.get_distance(diffx, diffy, boxx, boxy)
                if dist <= rcut:
                    atoms[i].neighbors.append(j)
                    atoms[j].neighbors.append(i)
                    atoms[i].diffxs.append(diffx)
                    atoms[j].diffxs.append(-diffx)
                    atoms[i].diffys.append(diffy)
                    atoms[j].diffys.append(-diffy)

        return atoms

    def write_traj(self, dump, i, box, atoms):
        """
        Single partcle dump function
        writes sno, x, y, vx, vy
        """
        #now write
        dump.write("ITEM: TIMESTEP\n")
        dump.write("&d\n"%i)
        dump.write("ITEM: NUMBER OF ATOMS\n")
        dump.write("%d\n" % len(atoms))
        dump.write("ITEM: BOX BOUNDS\n")
        dump.write("%f %f\n" % (box[0][0], box[0][1]))
        dump.write("%f %f\n" % (box[1][0], box[1][1]))
        dump.write("%f %f\n" % (-1, 1))
        title_str = "ITEM: ATOMS id type x y z vx vy vz\n"
        dump.write(title_str)
        for cc, atom in enumerate(atoms):
            atomline = ("%d %d %f %f %f %f %f %f\n")%(atom.id, atom.type, atom.x, atom.y, 0.0, atom.vx, atom.vy, 0.0)
            dump.write(atomline)
        dump.flush()

    def read_file(self, filename):
        """
        Read a dump file in column format
        Takes filename as input and returns atoms
        """
        atoms = []
        paramsread = False
        f = open(filename,'r')

        for count, line in enumerate(f):
            if not paramsread:
                #atom numer is at line 3
                if count == 3:
                    natoms = int(line.strip())
                #box dims in lines 5,6,7
                elif count == 5:
                    raw = line.strip().split()
                    boxx = [float(raw[0]), float(raw[1])]
                elif count == 6:
                    raw = line.strip().split()
                    boxy = [float(raw[0]), float(raw[1])]
                elif count == 7:
                    raw = line.strip().split()
                    boxz = [float(raw[0]), float(raw[1])]

                #header is here
                elif count == 8:
                    raw = line.strip().split()
                    headerdict = { raw[x]:x-2 for x in range(0, len(raw)) }
                    paramsread = True

            else:
                raw = line.strip().split()
                typ = int(raw[headerdict["type"]])
                x = float(raw[headerdict["x"]])
                y = float(raw[headerdict["y"]])
                vx = float(raw[headerdict["vx"]])
                vy = float(raw[headerdict["vy"]])

                atom = Atom(x, y, vx, vy)
                atom.type = typ

                atoms.append(atom)

        #close files
        f.close()
        return atoms, box
