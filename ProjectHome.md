## Python implementation of the fast radial basis function (RBF) interpolation for scientific applications. ##
Radial basis function (RBF) interpolation is a technique for representing a function starting with data on scattered points. Solving large RBF interpolation problems is notoriously difficult with basis functions of global support, due to the need to solve a linear system with a fully populated and badly conditioned matrix. Compact support bases result in sparse matrices, but at the cost of reduced approximation qualities.

In the present method, we keep the global basis functions, but localize their effect during the solution procedure, and add back the global effect of the bases via the iterations. The complexity of the present method is O(N) for N particles.

**March 2009** -- at this time, we release a stable alfa version of the Python code, and invite interested parties to email us if they would like to collaborate with us on further developments.

|We distribute this code under the MIT License, giving potential users the greatest freedom possible. We do, however, request fellow scientists that if they use our codes in research, they kindly include us in the acknowledgement of their papers.  We do not request gratuitous citations;  only cite our articles if you deem it warranted.|
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|

Visit the [group's webpage](http://barbagroup.bu.edu/) for more codes.

### Publications ###

> "Fast Radial Basis Function Interpolation with Gaussians by Localization and Iteration" by Claudio E. Torres and L. A. Barba

> _Journal of Computational Physics_, 2009 [DOI](http://dx.doi.org/10.1016/j.jcp.2009.03.007)

> Radial basis function (RBF) interpolation is a technique for representing a function starting with data on scattered points. In the vortex particle method for solving the Navier-Stokes equations, the representation of the ﬁeld of interest (the ﬂuid vorticity) as a sum of Gaussians is like an RBF approximation. In this application, there are two instances when one may need to solve an RBF interpolation problem: upon initialization of the vortex particles, and after spatial adaptation to replace a set of particles by another with a more regular distribution. Solving large RBF interpolation problems is notoriously difficult with basis functions of global support, due to the need to solve a linear system with a fully populated and badly conditioned matrix. In the vortex particle method one uses Gaussians of very small spread, which lead us to formulate a method of solution consisting of localization of the global problem, and improvement of solutions over iterations. It can be thought of as an opposite approach from the use of compact support basis functions to ease the solution procedure. Compact support bases result in sparse matrices, but at the cost of reduced approximation qualities. Instead, we keep the global basis functions, but localize their effect during the solution procedure, and add back the global effect of the bases via the iterations. Numerical experiments show that convergence is always obtained when the localized domains overlap via a buffer layer, and the particles overlap moderately. Algorithmic efficiency is excellent, achieving convergence to almost machine precision very rapidly with well-chosen parameters. We also discuss the complexity of the method and show that it can scale as O(N ) for N particles.
