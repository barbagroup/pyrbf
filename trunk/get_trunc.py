from pylab import *
def get_trunc(particle,cluster,ic):
    cluster.neighbor_trunc = int(ceil((cluster.nsigma_trunc-cluster.nsigma_box)/2/\
                          cluster.nsigma_box))
    cluster.trunc_length = cluster.nsigma_trunc*particle.sigma/2+1e-5

# loop through all clusters
    ista = cluster.ista[ic]
    iend = cluster.iend[ic]
    if ista <= iend:
        xc = cluster.xc[ic]
        yc = cluster.yc[ic]
        ix = cluster.ix[ic]
        iy = cluster.iy[ic]
        jx_min = max(0,ix-cluster.neighbor_trunc)
        jx_max = min(cluster.nx-1,ix+cluster.neighbor_trunc)
        jy_min = max(0,iy-cluster.neighbor_trunc)
        jy_max = min(cluster.ny-1,iy+cluster.neighbor_trunc)

# put all particles in the center box into the corresponding cell structure
        jsta = cluster.jsta[ic]
        jend = cluster.jend[ic]
        clusterx = []
        clustery = []
        clusterg = []
        clusterx.extend(particle.xj[jsta:jend+1])
        clustery.extend(particle.yj[jsta:jend+1])
        clusterg.extend(particle.gj[jsta:jend+1])

# loop through all neighbor
        for jx in range(jx_min,jx_max+1):
            for jy in range(jy_min,jy_max+1):
                if ix != jx or iy != jy:
                    jc = jx*cluster.ny+jy
                    jsta = cluster.jsta[jc]
                    jend = cluster.jend[jc]

# select from the particles in the neighbor boxes, the ones that belong in the trunc zone
                    xj = particle.xj[jsta:jend+1]
                    yj = particle.yj[jsta:jend+1]
                    gj = particle.gj[jsta:jend+1]
                    trunc = abs(xj-xc)<cluster.trunc_length
                    xj = extract(trunc,xj)
                    yj = extract(trunc,yj)
                    gj = extract(trunc,gj)
                    trunc = abs(yj-yc)<cluster.trunc_length
                    xj = extract(trunc,xj)
                    yj = extract(trunc,yj)
                    gj = extract(trunc,gj)
                
# add all particles in the neighbor boxes into the corresponding cell structure
                    clusterx.extend(xj)
                    clustery.extend(yj)
                    clusterg.extend(gj)
            cluster.xjt=clusterx
            cluster.yjt=clustery
            cluster.gjt=clusterg
            cluster.nptruncj = len(clusterx)
    else:
        cluster.xjt=[]
        cluster.yjt=[]
        cluster.gjt=[]
        cluster.nptruncj = 0

    return particle,cluster
