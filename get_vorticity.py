from pylab import *
from get_trunc import *
def get_vorticity(particle,cluster):
    particle.wi = arange(particle.ni)*1.
    for ic in range(cluster.n):
        particle,cluster = get_trunc(particle,cluster,ic)
        ista = cluster.ista[ic]
        iend = cluster.iend[ic]
        for ip in range(ista,iend+1):
            dx = particle.xi[ip]-cluster.xjt
            dy = particle.yi[ip]-cluster.yjt
            particle.wi[ip] = dot(cluster.gjt,exp(-(dx**2+dy**2)/\
                (2*particle.sigma**2))/(2*pi*particle.sigma**2))
    return particle,cluster
