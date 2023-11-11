import numpy as np
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
RESTART_FILE = 0
RC_FNAME = "pgrid4py.dat" ## coord file to read PREFIX
D_FNAME =  "pygrid_tot.dat"    ## file to write PREFIX

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
print("Making clean initial conditions")

### check if restart

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
