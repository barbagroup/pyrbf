from pylab import *
def get_buffer(particle,cluster,ic):
    cluster.neighbor_buffer = int(ceil((cluster.nsigma_buffer-cluster.nsigma_box)/2/\
                          cluster.nsigma_box))
    cluster.buffer_length = cluster.nsigma_buffer*particle.sigma/2+1e-5

# loop through all i clusters
    ista = cluster.ista[ic]
    iend = cluster.iend[ic]
    if ista <= iend:
        xc = cluster.xc[ic]
        yc = cluster.yc[ic]
        ix = cluster.ix[ic]
        iy = cluster.iy[ic]
        jx_min = max(0,ix-cluster.neighbor_buffer)
        jx_max = min(cluster.nx-1,ix+cluster.neighbor_buffer)
        jy_min = max(0,iy-cluster.neighbor_buffer)
        jy_max = min(cluster.ny-1,iy+cluster.neighbor_buffer)
    
# put all particles in the center box into the corresponding cell structure
        clusterx = []
        clustery = []
        clusterg = []
        clustere = []
        clusterw = []
        clusterx.extend(particle.xi[ista:iend+1])
        clustery.extend(particle.yi[ista:iend+1])
        clusterg.extend(particle.gj[ista:iend+1])
        clustere.extend(particle.ei[ista:iend+1])
        clusterw.extend(particle.wi[ista:iend+1])
    
# loop through all neighbor
        for jx in range(jx_min,jx_max+1):
            for jy in range(jy_min,jy_max+1):
                if ix != jx or iy != jy:
                    jc = jx*cluster.ny+jy
                    jsta = cluster.ista[jc]
                    jend = cluster.iend[jc]
                
# select from the particles in the neighbor boxes, the ones that belong in the buffer zone
                    xi = particle.xi[jsta:jend+1]
                    yi = particle.yi[jsta:jend+1]
                    gi = particle.gj[jsta:jend+1]
                    ei = particle.ei[jsta:jend+1]
                    wi = particle.wi[jsta:jend+1]
                    buffer = abs(xi-xc)<cluster.buffer_length
                    xi = extract(buffer,xi)
                    yi = extract(buffer,yi)
                    gi = extract(buffer,gi)
                    ei = extract(buffer,ei)
                    wi = extract(buffer,wi)
                    buffer = abs(yi-yc)<cluster.buffer_length
                    xi = extract(buffer,xi)
                    yi = extract(buffer,yi)
                    gi = extract(buffer,gi)
                    ei = extract(buffer,ei)
                    wi = extract(buffer,wi)
                
# add all particles in the neighbor boxes into the corresponding cell structure
                    clusterx.extend(xi)
                    clustery.extend(yi)
                    clusterg.extend(gi)
                    clustere.extend(ei)
                    clusterw.extend(wi)
            cluster.xib=clusterx
            cluster.yib=clustery
            cluster.gib=clusterg
            cluster.eib=clustere
            cluster.wib=clusterw
            cluster.npbufferi = len(clusterx)
    else:
        cluster.xib=[]
        cluster.yib=[]
        cluster.gib=[]
        cluster.eib=[]
        cluster.wib=[]
        cluster.npbufferi = 0
            
    return particle,cluster