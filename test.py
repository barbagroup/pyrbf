# Fast radial basis function interpolation for 2-D Lamb-Oseen vortex
# Rio Yokota 2009/03/16
#
# Basic convention for variable names
# parameter.* : physics oriented parameters
# particle.*  : particle oriented variables
# *   : cluster oriented variables
# x           : x coordinate
# y           : y coordinate
# g           : vortex strength \gamma
# e           : exact vorticity \omega
# w           : calculated vorticity \omega_h

from pylab import *
from rbf_solver2 import *

epsf = 1e-8
epsd = 1e-20

# physical parameters
vis = 0.1; t = 1
xmin = -1; xmax = 1
ymin = -1; ymax = 1
solver = 0; wrapper = 1;

# particle parameters
sigma = 0.05; overlap = 1.0
h = overlap*sigma

# generate particles
x = arange(xmin,xmax+epsf,h)
y = arange(ymin,ymax+epsf,h)
nx = x.shape[0]
ny = y.shape[0]
x,y = meshgrid(x,y)
x = x.flatten(1); y = y.flatten(1)
n = nx*ny
    
# calculate initial vortex strength and exact vorticity field on particle
e = exp(-(x**2+y**2)/\
    (4*vis*t))/(pi*4*vis*t)
g = e*h**2
w = e

# solve rbf using c wrapper
if wrapper == 1:
    rbf_solver(x,y,g,e,sigma)
