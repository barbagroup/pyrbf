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
from fmm import *

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
ni = nx*ny
nj = ni
    
# calculate initial vortex strength and exact vorticity field on particle
e = exp(-(x**2+y**2)/(4*vis*t))/(pi*4*vis*t)
g = e*h**2
w = e

np = 2e6
xi = zeros(np)
yi = zeros(np)
ui = zeros(np)
vi = zeros(np)
xj = zeros(np)
yj = zeros(np)
gj = zeros(np)
sj = zeros(np)
for i in range(ni):
    xi[i] = x[i]
    yi[i] = y[i]
    xj[i] = x[i]
    yj[i] = y[i]
    gj[i] = g[i]
    sj[i] = sigma

bs2(1,ni,xi,yi,ui,vi,1,nj,xj,yj,gj,sj,10,5)
for i in range(ni):
    uij = 0.0 
    vij = 0.0 
    for j in range(nj):
        dxij = xi[i]-xj[j]
        dyij = yi[i]-yj[j]
        rij = dxij**2+dyij**2+epsd
        sij = sj[j]**2
        cutoff = 1/rij*(1-exp(-rij/2/sij))
        uij = uij+0.5/pi*gj[j]*dyij*cutoff
        vij = vij-0.5/pi*gj[j]*dxij*cutoff
    print i,ui[i],uij,vi[i],vij
