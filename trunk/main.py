# Fast radial basis function interpolation for 2-D Lamb-Oseen vortex
# Rio Yokota 2009/03/16
#
# Basic convention for variable names
# parameter.* : physics oriented parameters
# particle.*  : particle oriented variables
# cluster.*   : cluster oriented variables
# x           : x coordinate
# y           : y coordinate
# g           : vortex strength \gamma
# e           : exact vorticity \omega
# w           : calculated vorticity \omega_h

from pylab import *
from resource import *
from get_cluster import *
from get_vorticity import *
from get_strength import *
from rbf_solver2 import rbf_solver

def cpu_time():
    return (getrusage(RUSAGE_SELF).ru_utime+
            getrusage(RUSAGE_SELF).ru_stime)
tic = cpu_time()
t = arange(10)*0.

# define classes
class PARAMETER:
    pass
class PARTICLE:
    pass
class CLUSTER:
    pass
class GRID:
    pass
class SOLVER:
    pass
parameter = PARAMETER()
particle = PARTICLE()
cluster = CLUSTER()
grid = GRID()
solver = SOLVER()
epsf = 1e-8
epsd = 1e-20

# physical parameters
parameter.vis = 0.1; parameter.t = 1
solver.method = 0; parameter.wrapper = 0

# particle parameters
particle.sigma = 0.05; particle.overlap = 1.0
particle.h = particle.overlap*particle.sigma
particle.xmin = 0; particle.xmax = 1
particle.ymin = 0; particle.ymax = 1

# cluster parameters
cluster.nsigma_box = 8.
cluster.nsigma_buffer = 12.
cluster.nsigma_trunc = 24.

# generate particles
x = arange(particle.xmin,particle.xmax+epsf,particle.h)
y = arange(particle.ymin,particle.ymax+epsf,particle.h)
grid.nx = x.shape[0]
grid.ny = y.shape[0]
grid.x,grid.y = meshgrid(x,y)
particle.xi = grid.x.flatten(1); particle.yi = grid.y.flatten(1)
particle.xj=particle.xi.copy(); particle.yj=particle.yi.copy()
particle.ni = grid.nx*grid.ny
particle.nj=particle.ni
    
# calculate initial vortex strength and exact vorticity field on particle
#particle.ei = exp(-(particle.xi**2+particle.yi**2)/\
#    (4*parameter.vis*parameter.t))/(pi*4*parameter.vis*parameter.t))
particle.ei = 0.75*exp(-((9*particle.xi-2)**2+(9*particle.yi-2)**2)/4)\
    +0.75*exp(-((9*particle.xi+1)**2)/49-((9*particle.yi+1)**2)/10)\
    +0.5*exp(-((9*particle.xi-7)**2+(9*particle.yi-3)**2)/4)\
    -0.2*exp(-(9*particle.xi-4)**2-(9*particle.yi-7)**2)
particle.wi = particle.ei
#particle.gj = exp(-(particle.xj**2+particle.yj**2)/\
#    (4*parameter.vis*parameter.t))/(pi*4*parameter.vis*parameter.t)*particle.h**2
particle.gj = (0.75*exp(-((9*particle.xj-2)**2+(9*particle.yj-2)**2)/4)\
    +0.75*exp(-((9*particle.xj+1)**2)/49-((9*particle.yj+1)**2)/10)\
    +0.5*exp(-((9*particle.xj-7)**2+(9*particle.yj-3)**2)/4)\
    -0.2*exp(-(9*particle.xj-4)**2-(9*particle.yj-7)**2))*particle.h**2

toc = tic
tic = cpu_time()
t[0] = tic-toc

# solve rbf using c wrapper
if parameter.wrapper == 1:
    rbf_solver(particle.xi,particle.yi,particle.gj,particle.ei,
        particle.sigma,parameter.vis)

toc = tic
tic = cpu_time()

# generate clusters
particle,cluster = get_cluster(particle,cluster)

# set up reference data
if parameter.wrapper == 1:
    particle,cluster = get_trunc(particle,cluster)
    particle,cluster = get_vorticity(particle,cluster)
    particle.gj = particle.ei*particle.h**2
    particle.ei = particle.wi
    particle.wi = particle.gj/particle.h**2

toc = tic
tic = cpu_time()
t[0] += tic-toc

# RBF interpolation
it = -1; iconv = 0; err = []; grid.r = []
while iconv < 5:
    it = it+1

# estimate vorticity field on particle from vortex strength
    particle,cluster = get_vorticity(particle,cluster)

    toc = tic
    tic = cpu_time()
    t[1] += tic-toc
    
# solve the system of equations to calculate the vortex strength
    particle,cluster,solver = get_strength(particle,cluster,solver,it)

    toc = tic
    tic = cpu_time()
    t[2] += tic-toc

# calculate the L2 norm error
    err.append(log10(std(solver.r)/std(particle.ei)))
    print 'iteration : ',it,' error : ',err[it]
    if err[it] < -14: iconv = iconv+1
    
# plot the spatial variation of error
    particle.ri[particle.isort] = particle.ri.copy()
    grid.r.append(particle.ri.reshape(grid.ny,-1))
    imshow(log10(abs(grid.r[it])/max(abs(particle.ei))+epsd))
    if it == 0: colorbar()
    clim(-20, 0); draw()

    toc = tic
    tic = cpu_time()
    t[0] += tic-toc

for i in range(3): t[9] += t[i]
print 'matvec : ',t[1]
print 'solver : ',t[2]
print 'other  : ',t[0]
print '------------------'
print 'total  : ',t[9]
