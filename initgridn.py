import numpy as np
import sys

prank = sys.argv[1] ##need to supply rank of process

grid_fname = "pygrid_tot.dat"
d_fname = f"pygrid_{int(prank):03d}.dat"    ## file to write PREFIX
ngrid_fname = f"pgrid4py_{int(prank):03d}.dat"
####################################################################################################
###
### MIN_PRS:   Sets background to some minimum pressure, specified as a fraction of the maximum
### MIN_RHO:   Sets background to some minimum rho, specified as a fraction of the maximum
### FORCE_GEQ: Sets the background to the density minimum and picks P to force radial equilibrium
###
####################################################################################################

# convenient definitions
prims = ["rho", "prs", "vx1", "vx2", "vx3"]
X1, X2, RHO, PRS, VX1, VX2, VX3 = range(0,7)

## need to read the total grid
gfile = open(grid_fname, "r")
rcs   = np.fromstring(gfile.readline(), sep="\t")
thcs  = np.fromstring(gfile.readline(), sep="\t")
NRCS  = len(rcs)
NTHCS = len(thcs)
cgrid = np.zeros((NRCS, NTHCS, 7))

## reconstruct total cgrid

for j, th in enumerate(thcs):
    for i,r in enumerate(rcs):
        cgrid[i,j,0] = r
        cgrid[i,j,1] = th
        cgrid[i,j, 2:] = np.fromstring(gfile.readline(), sep="\t")

gfile.close()
print(f"Total grid: RCOORDS({NRCS})  = [{rcs[0]},{rcs[-1]}]",
                   f"THCOORDS({NTHCS}) = [{thcs[0]},{thcs[-1]}]")

## load subgrid
ngfile = open(ngrid_fname, "r")
nrcs = np.fromstring(ngfile.readline(),sep="\t")
nthcs = np.fromstring(ngfile.readline(),sep="\t")
nNRCS = len(nrcs)
nNTHCS = len(nthcs)
ngfile.close()
print(f"Subgrid: RCOORDS({nNRCS})  = [{nrcs[0]},{nrcs[-1]}]",
                   f"THCOORDS({nNTHCS}) = [{nthcs[0]},{nthcs[-1]}]")

# defcgrid = np.zeros((NRCS, NTHCS, 7))
####################################################################################################
###
### defines that function that gives a (the closest) grid coordinate give a spatial point 
###
####################################################################################################

def c2index(r,th):
    rdiff = np.abs(rcs-r)
    thdiff = np.abs(thcs - th)

    if np.min(rdiff) > 1e-3:
        print(f"r coordinate imprecise {np.min(rdiff):5.5e}")
    if np.min(thdiff) > 1e-3:
        print(f"th coordinate imprecise {np.min(thdiff):5.5e}")
    
    rindex = np.argmin(rdiff)
    thindex= np.argmin(thdiff)
    return (rindex, thindex)

def write_ncgrid():
    ncwrite = open(d_fname ,"w")
    for j,th in enumerate(nthcs):
        for i,r in enumerate(nrcs):
           cgi, cgj = c2index(r,th)
           rho,prs,vx1,vx2,vx3 = cgrid[cgi,cgj,2:]
           ncwrite.write(f"{rho:5.5e}\t{prs:5.5e}\t{vx1:5.5e}\t{vx2:5.5e}\t{vx3:5.5e}\n")
    ncwrite.close()
    print(f"Wrote file {d_fname}")

####################################################################################################
##
## write files
##
####################################################################################################
write_ncgrid()

