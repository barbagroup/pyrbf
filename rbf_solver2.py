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
from get_buffer import *
from get_trunc import *
from get_vorticity import *
from get_strength import *

def cpu_time():
    return (getrusage(RUSAGE_SELF).ru_utime+getrusage(RUSAGE_SELF).ru_stime)
    
    
def rbf_solver(x,y,g,e,sigma):

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
    particle.ni = x.shape[0]
    particle.nj = particle.ni
    particle.sigma = sigma
    particle.xmin = x.min()
    particle.xmax = x.max()
    particle.ymin = y.min()
    particle.ymax = y.max()
    solver.method = 0

# cluster parameters
    cluster.nsigma_box = 8.
    cluster.nsigma_buffer = 12.
    cluster.nsigma_trunc = 24.

# initialize particles
    particle.xi = x
    particle.yi = y
    particle.ei = e
    particle.wi = e
    particle.xj = x
    particle.yj = y
    particle.gj = g
    
# generate clusters
    particle,cluster = get_cluster(particle,cluster)

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
        err.append(log10(std(particle.ri)/std(particle.ei)))
        print 'iteration : ',it,' error : ',err[it]
        if err[it] < -5: iconv = iconv+1

        toc = tic
        tic = cpu_time()
        t[0] += tic-toc

    for i in range(3): t[9] += t[i]
    print 'matvec : ',t[1]
    print 'solver : ',t[2]
    print 'other  : ',t[0]
    print '------------------'
    print 'total  : ',t[9]
