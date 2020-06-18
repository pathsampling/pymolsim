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

class SingleParticle:
    """
    Class contains single particle potentials and methods to write/read trajectories
    """
    def __init__(self):
        pass

    def simple_2d(self, x, y, **kwargs):
        """
        Potential for simple2d with modifyable parameters.

        """
        g1 = kwargs.get('g1', 1.00)
        g2 = kwargs.get('g2', -10.00)
        g3 = kwargs.get('g3', -10.00)
        a1 = kwargs.get('a1', -30.00)
        a2 = kwargs.get('a2', -3.00)
        b1 = kwargs.get('b1', -30.00)
        b2 = kwargs.get('b2', -3.00)
        x0 = kwargs.get('x0', 0.2)
        y0 = kwargs.get('y0', 0.0)

        t1     = g1*(x*x + y*y)**2
        dt1_dx = g1*2*(x*x + y*y)*2*x
        dt1_dy = g1*2*(x*x + y*y)*2*y

        t2     = g2*np.exp( a1*(x-x0)**2 + a2*(y-y0)**2)
        dt2_dx = t2*a1*2*(x-x0)
        dt2_dy = t2*a2*2*(y-y0)

        t3     = g3*np.exp( b1*(x+x0)**2 + b2*(y+y0)**2)
        dt3_dx = t3*b1*2*(x+x0)
        dt3_dy = t3*b2*2*(y+y0)

        potential = t1 + t2 + t3
        force_x = -(dt1_dx + dt2_dx + dt3_dx)
        force_y = -(dt1_dy + dt2_dy + dt3_dy)

        return potential, force_x, force_y

    def write_traj(self, dump, dumpdata):
        """
        Single partcle dump function
        writes sno, x, y, vx, vy
        """
        dump.write("%d %f %f %f %f\n"%(dumpdata[0],dumpdata[1],dumpdata[2],dumpdata[3],dumpdata[4]))
        dump.flush()

    def read_file(self, filename):
        """
        Read a dump file in column format
        Takes filename as input and returns atoms
        """
        atoms = []
        idd, x, y, vx, vy = np.loadtxt(filename, unpack=True, dtype=float)
        at = Atom(x, y, vx, vy)
        atoms.append(at)

        return atoms
