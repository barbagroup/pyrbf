from pylab import *
from get_buffer import *
# domain decomposition with direct solver
def DD(particle,cluster,solver,it):
    solver.r = particle.ei-particle.wi
    for ic in range(cluster.n):
        particle,cluster = get_buffer(particle,cluster,ic)
        ista = cluster.ista[ic]
        iend = cluster.iend[ic]
        if ista <= iend:
            [xi,xj] = meshgrid(cluster.xib,cluster.xib)
            dx = xi-xj
            [yi,yj] = meshgrid(cluster.yib,cluster.yib)
            dy = yi-yj
            A = exp(-(dx**2+dy**2)/(2*particle.sigma**2))/(2*pi*particle.sigma**2)
            b = array(cluster.eib)-array(cluster.wib)+dot(array(cluster.gib),A)
            x = solve(A,b)
            particle.gj[ista:iend+1] = x[0:iend-ista+1]
    return particle,cluster,solver
# CG
def CG(particle,cluster,solver,it):
    if it == 0:
        solver.r = particle.ei.copy()
        solver.x = solver.r*0.
        solver.p = solver.r.copy()
        solver.rho = dot(solver.r,solver.r)
        particle.gj = solver.p.copy()
    solver.q = particle.wi.copy()
    solver.alpha = solver.rho/dot(solver.p,solver.q)
    solver.x += solver.alpha*solver.p
    solver.r -= solver.alpha*solver.q
    solver.rho_old = solver.rho
    solver.rho = dot(solver.r,solver.r)
    solver.beta = solver.rho/solver.rho_old
    solver.p = solver.r+solver.beta*solver.p
    particle.gj = solver.p.copy()
    return particle,cluster,solver

# BICG
def BICG(particle,cluster,solver,it):
    if it == 0:
        solver.r = particle.ei.copy()
        solver.rt = solver.r.copy()
        solver.x = solver.r*0.
        solver.p = solver.r.copy()
        solver.pt = solver.rt.copy()
        solver.rho = dot(solver.r,solver.rt)
        particle.gj = solver.p.copy()
    if it%2 == 0:
        solver.q = particle.wi.copy()
        solver.alpha = solver.rho/dot(solver.pt,solver.q)
        solver.x += solver.alpha*solver.p
        solver.r -= solver.alpha*solver.q
        particle.gj = solver.pt.copy()
    else:
        solver.qt = particle.wi.copy()
        solver.rt -= solver.alpha*solver.qt
        solver.rho_old = solver.rho
        solver.rho = dot(solver.r,solver.rt)
        solver.beta = solver.rho/solver.rho_old
        solver.p = solver.r+solver.beta*solver.p
        solver.pt = solver.rt+solver.beta*solver.pt
        particle.gj = solver.p.copy()
    return particle,cluster,solver

#BICGSTAB
def BICGSTAB(particle,cluster,solver,it):
    if it == 0:
        solver.r = particle.ei.copy()
        solver.rt = solver.r.copy()
        solver.x = solver.r*0.
        solver.p = solver.r.copy()
        solver.rho = dot(solver.r,solver.rt)
        particle.gj = solver.p.copy()
    if it%2 == 0:
        solver.v = particle.wi.copy()
        solver.alpha = solver.rho/dot(solver.rt,solver.v)
        solver.s = solver.r-solver.alpha*solver.v
        if std(solver.s)<1e-2:
            solver.x += solver.alpha*solver.p
        particle.gj = solver.s.copy()
    else:
        solver.t = particle.wi.copy()
        solver.omega = dot(solver.t,solver.s)/dot(solver.t,solver.t)
        solver.x += solver.alpha*solver.p+solver.omega*solver.s
        solver.r = solver.s-solver.omega*solver.t
        solver.rho_old = solver.rho
        solver.rho = dot(solver.r,solver.rt)
        solver.beta = solver.rho/solver.rho_old*solver.alpha/solver.omega
        solver.p = solver.r+solver.beta*(solver.p-solver.omega*solver.v)
        particle.gj = solver.p.copy()
    return particle,cluster,solver
def get_strength(particle,cluster,solver,it):
    if solver.method == 0:
        particle,cluster,solver = DD(particle,cluster,solver,it)
    elif solver.method == 1:
        particle,cluster,solver = CG(particle,cluster,solver,it)
    elif solver.method == 2:
        particle,cluster,solver = BICG(particle,cluster,solver,it)
    elif solver.method == 3:
        particle,cluster,solver = BICGSTAB(particle,cluster,solver,it)
    particle.ri=solver.r.copy()
    return particle,cluster,solver