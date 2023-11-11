import numpy as np
import sys, nsim
from scipy.integrate import dblquad

NGAM = 3.       # polytropic index
MBH  = 2.e5     # mass of central object
MDISK= 2.e3     # total mass of disk
RDISK= 5        # radius of disk in units that we use
DIST = 1.04     # distortion parameter
RHOCRIT = 1.e-4      # minimum density as a fraction of maximum density
PRSMIN = 1.e-8       # minimum value of pressure (doesn't really matter)
GAMMA = 1 + 1/NGAM
GRAV_CONST = 6.67e-2 
SPECAM     = np.sqrt(GRAV_CONST*MBH*RDISK)

RESTART_DIR = "restart"
RESTART_FILE = 1
RESTART_INDEX = np.inf
RESTART_MODVEL = True ## reset non azimuthal velocities
RESTART_INDEX  = 3
RC_FNAME = "pgrid4py.dat" ## coord file to read PREFIX
D_FNAME =  "pygrid_tot.dat"    ## file to write PREFIX
NGHOST = 3 ## number of ghost zones
#######################################################################
## need to do some stuff for parallelization 
#####################################################################

f = open(RC_FNAME, "r")
rcs = np.sort(np.fromstring(f.readline(), sep="\t"))
thcs= np.sort(np.fromstring(f.readline(), sep="\t"))
f.close()
NTHCS= len(thcs)
NRCS = len(rcs)
print(f"Total grid: RCOORDS({NRCS})  = [{rcs[0]},{rcs[-1]}]",
                   f"THCOORDS({NTHCS}) = [{thcs[0]},{thcs[-1]}]")

# convenient definitions
prims = ["rho", "prs", "vx1", "vx2", "vx3"]
X1, X2, RHO, PRS, VX1, VX2, VX3 = range(0,7)

# defining the grid
cgrid = np.zeros((NRCS, NTHCS, 7))

### check if restart
if RESTART_FILE:
    print("Restarting file")
    sim = nsim.Simulation(RESTART_DIR)
    vx1 = sim.P["vx1"][RESTART_INDEX]
    vx2 = sim.P["vx1"][RESTART_INDEX]
    vx3 = sim.P["vx3"][RESTART_INDEX]
    rho = sim.P["rho"][RESTART_INDEX]
    prs = sim.P["prs"][RESTART_INDEX]
    
    ## populate center of cgrid
    for j, th in enumerate(thcs[NGHOST:NTHCS-NGHOST]):
        for i, r in enumerate(rcs[NGHOST:NRCS-NGHOST]):
            J = j + NGHOST
            I = i + NGHOST
            cgrid[I,J,RHO] = rho[i,j]
            cgrid[I,J,PRS] = prs[i,j]
            cgrid[I,J,VX3] = vx3[i,j]
            cgrid[I,J,VX1] = 0
            cgrid[I,J,VX2] = 0
        ## do the radial ghost zones
        for i in range(0, NGHOST):
            ## inner
            cgrid[i,J,RHO] = cgrid[NGHOST,J,RHO]
            cgrid[i,J,PRS] = cgrid[NGHOST,J,PRS]
            cgrid[i,J,VX1] = cgrid[NGHOST,J,VX1]
            cgrid[i,J,VX2] = cgrid[NGHOST,J,VX2]
            cgrid[i,J,VX3] = cgrid[NGHOST,J,VX3]

            cgrid[NRCS-i-1,J,RHO] = cgrid[NRCS-NGHOST-1,J,RHO]
            cgrid[NRCS-i-1,J,PRS] = cgrid[NRCS-NGHOST-1,J,PRS]
            cgrid[NRCS-i-1,J,VX1] = cgrid[NRCS-NGHOST-1,J,VX1]
            cgrid[NRCS-i-1,J,VX2] = cgrid[NRCS-NGHOST-1,J,VX2]
            cgrid[NRCS-i-1,J,VX3] = cgrid[NRCS-NGHOST-1,J,VX3]

    ### now do theta ones
    for i, r in enumerate(rcs):
        for j in range(0, NGHOST):
            cgrid[i,j,RHO] = cgrid[i,NGHOST,RHO]
            cgrid[i,j,PRS] = cgrid[i,NGHOST,PRS]
            cgrid[i,j,VX1] = cgrid[i,NGHOST,VX1]
            cgrid[i,j,VX2] = cgrid[i,NGHOST,VX2]
            cgrid[i,j,VX3] = cgrid[i,NGHOST,VX3]

            cgrid[i,NTHCS-1-j,RHO] = cgrid[i,NTHCS-NGHOST-1,RHO]
            cgrid[i,NTHCS-1-j,PRS] = cgrid[i,NTHCS-NGHOST-1,PRS]
            cgrid[i,NTHCS-1-j,VX1] = cgrid[i,NTHCS-NGHOST-1,VX1]
            cgrid[i,NTHCS-1-j,VX2] = cgrid[i,NTHCS-NGHOST-1,VX2]
            cgrid[i,NTHCS-1-j,VX3] = cgrid[i,NTHCS-NGHOST-1,VX3]
else:
    print("making clean initial conditions")
    ## Mass integral stuff
    def Idn_integrand(r, th):
        return (r**2*np.sin(th)) * (1/r-0.5*(1/r/np.sin(th))**2 - 0.5/DIST) ** NGAM

    ## defining bounds of integration
    def r_lower(th):
        return (np.sin(th)**2*(1+np.sqrt(1-1/DIST/np.sin(th)**2)))**-1

    def r_upper(th):
        return (np.sin(th)**2*(1-np.sqrt(1-1/DIST/np.sin(th)**2)))**-1

    th_upper = np.arcsin(DIST**-0.5)

    Idn = 2*dblquad(Idn_integrand, th_upper, np.pi/2, r_lower, r_upper)[0]

    alpha = 2*np.pi*(GRAV_CONST*MBH/(NGAM+1))**NGAM*Idn
    entA  = (MDISK/(alpha * RDISK**(3-NGAM)))**(-1/NGAM)
    RHOMIN     = RHOCRIT*(GRAV_CONST*MBH/2/(NGAM+1)/entA/RDISK*(1-1/DIST))**NGAM

    ####################################################################################################
    ###
    ### Sets up the mass profile (doesn't do anything about the ambient medium)
    ###
    ####################################################################################################

    def rho0(r, th):
        const = (GRAV_CONST * MBH/(NGAM+1)/entA/RDISK)**NGAM
        cdstuff = (RDISK/r-0.5*(RDISK/(r*np.sin(th)))**2-0.5/DIST)**NGAM
        return const*cdstuff

    def phi(r):
        return -GRAV_CONST*MBH/r

    ### initializes grid (no pressure)
    for i, r in enumerate(rcs):
        for j, th in enumerate(thcs):
            rho = rho0(r, th)
            cgrid[i,j,X1] = r
            cgrid[i,j,X2] = th
            if rho < RHOMIN:
                cgrid[i,j,RHO] = RHOMIN
                continue
            cgrid[i,j,RHO] = rho
            cgrid[i,j,VX3] = SPECAM/( r*np.sin(th) )

    ### Initializes the pressure
    for j, th in enumerate(thcs):
        for i,r in enumerate(rcs[::-1]):
            I=np.shape(cgrid)[0]-i-1
            if i == 0:
                cgrid[I,j,PRS] = 0
                continue
            dr = cgrid[I+1,j,X1] - cgrid[I,j,X1]
            dP = cgrid[I,j,RHO]*dr/cgrid[I,j,X1]*(cgrid[I,j,VX3]**2+phi(cgrid[I,j,X1]))
            cgrid[I,j,PRS] = cgrid[I+1,j,PRS] - dP
            cgrid[:,:,PRS] += 2*PRSMIN

            
    ####################################################################################################
    ###
    ### Writes to file
    ###
    ####################################################################################################


## input PCOORD for process -> output cgrid for pcoords
def write_cgrid():
    f   = open(D_FNAME, "w")
## print r and theta coordiantes to the file
    for i,r in enumerate(rcs):
        f.write(f"{r}\t")
        
    f.write("\n")
    for j, th in enumerate(thcs):
        f.write(f"{th}\t")
    f.write("\n")
    
    for j,th in enumerate(thcs):
        for i,r in enumerate(rcs):
           rho,prs,vx1,vx2,vx3 = cgrid[i,j,2:]
           f.write(f"{rho:5.5e}\t{prs:5.5e}\t{vx1:5.5e}\t{vx2:5.5e}\t{vx3:5.5e}\n")
    f.close()
    print(f"Wrote file {D_FNAME}")

####################################################################################################
##
## write files
##
####################################################################################################

write_cgrid()
