from pylab import *
def get_cluster(particle,cluster):
    epsf = 1e-8
    # calculate cluster size
    cluster.xmin = particle.xmin-epsf
    cluster.xmax = particle.xmax+epsf
    cluster.ymin = particle.ymin-epsf
    cluster.ymax = particle.ymax+epsf
    cluster.box_length = cluster.nsigma_box*particle.sigma+epsf
    
    # calculate number of clusters in each direction
    cluster.nx = int(ceil((cluster.xmax-cluster.xmin)/cluster.box_length))
    cluster.ny = int(ceil((cluster.ymax-cluster.ymin)/cluster.box_length))
    cluster.n = cluster.nx*cluster.ny

    # calculate the x,y index and coordinates of the center
    cluster.ix = arange(cluster.n)
    cluster.iy = arange(cluster.n)
    cluster.xc = arange(cluster.n)*1.
    cluster.yc = arange(cluster.n)*1.
    cluster.ista = arange(cluster.n)
    cluster.iend = arange(cluster.n)
    cluster.jsta = arange(cluster.n)
    cluster.jend = arange(cluster.n)
    ic = -1
    for ix in range(cluster.nx):
        for iy in range(cluster.ny):
            ic = ic+1
            cluster.ix[ic] = ix
            cluster.iy[ic] = iy
            cluster.xc[ic] = cluster.xmin+(ix+0.5)*cluster.box_length
            cluster.yc[ic] = cluster.ymin+(iy+0.5)*cluster.box_length
            cluster.ista[ic] = 0
            cluster.iend[ic] = -1
            cluster.jsta[ic] = 0
            cluster.jend[ic] = -1
    
    # assign cluster number to particles
    cluster.i = arange(particle.ni)
    cluster.j = arange(particle.nj)
    for ip in range(particle.ni):
        ix_cluster = int(floor((particle.xi[ip]-cluster.xmin)/cluster.box_length))
        iy_cluster = int(floor((particle.yi[ip]-cluster.ymin)/cluster.box_length))
        cluster.i[ip] = ix_cluster*cluster.ny+iy_cluster

    for ip in range(particle.nj):
        ix_cluster = int(floor((particle.xj[ip]-cluster.xmin)/cluster.box_length))
        iy_cluster = int(floor((particle.yj[ip]-cluster.ymin)/cluster.box_length))
        cluster.j[ip] = ix_cluster*cluster.ny+iy_cluster
        
    # sort particles according to cluster number
    particle.isort = cluster.i.argsort();
    cluster.i = cluster.i[particle.isort]
    particle.xi = particle.xi[particle.isort]
    particle.yi = particle.yi[particle.isort]
    particle.ei = particle.ei[particle.isort]
    particle.wi = particle.wi[particle.isort]
    particle.jsort = cluster.j.argsort();
    cluster.j = cluster.j[particle.jsort]
    particle.xj = particle.xj[particle.jsort]
    particle.yj = particle.yj[particle.jsort]
    particle.gj = particle.gj[particle.jsort]
    
    # calculate start and end index of the particles in each cluster
    ic = -1
    id = -1
    for ip in range(particle.ni):
        if id != cluster.i[ip]:
            ic += 1
            id = cluster.i[ip]
            cluster.ista[id] = ip
            if ic >= 1: cluster.iend[io] = ip-1
            io = id
    cluster.iend[id] = particle.ni-1

    ic = -1
    id = -1
    for ip in range(particle.nj):
        if id != cluster.j[ip]:
            ic += 1
            id = cluster.j[ip]
            cluster.jsta[id] = ip
            if ic >= 1: cluster.jend[io] = ip-1
            io = id
    cluster.jend[id] = particle.nj-1
    
    return particle,cluster