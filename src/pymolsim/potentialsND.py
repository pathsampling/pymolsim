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
from pathsampling.potentials1D import *
from pathsampling.potentials2D import *
from pathsampling.potentialsmd import *
