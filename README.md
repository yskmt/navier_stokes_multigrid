3D Navier Stokes Equation Solver                                                                                                                                        
================================

Solving 3D incompressible Navier Stokes equation using finite difference method with uniform grid in parallel. The incompressibility is implemented using pressure-corr\
ection scheme and linear (Poisson) solver is implemented using multigrid v-cycle.

The code is pretty much based on this [MATLAB implementation](http://math.mit.edu/cse/codes/mit18086_navierstokes.m). The documentation is [here](http://math.mit.edu/c\
se/codes/mit18086_navierstokes.pdf).

Problem Statement
-------------------------

Incompressible Navier Stokes Equation:

![momentum](http://upload.wikimedia.org/math/4/8/c/48c88ec1a44dce97a23ceff09ee668b2.png "momentum")

Incompressibility condition:

![mass](http://upload.wikimedia.org/math/1/6/9/169892c54316eb6d350f5118bff5c213.png "mass")

3D cubic domain [0,1]^3. 

Dicirhlet boundary conditions for veocities and Neumann boundary condition for pressure are implemented. 
The velocity on the boundary can be specified in double bcs[3][6] array (first dimension specifies x- y- z- velocity,
second dimsension specifies the face of the cube).

Implementation
-------------------------

1. The cubic domain is discretized by unform staggered grid.
2. The partial differential equations are discretized using finite difference method.
3. Parallelization is implemented by OpenMP.
4. Pressure-correction scheme is used to enforce incompressibility.

Usage
-------------------------

./multigrid [# of threads] [max level of v-cycle] [x-grid size] [y-grid size] [z-grid size]


