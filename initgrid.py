import numpy as np
from scipy.integrate import dblquad
import glob
import sys
import plutotools as pt

## for restarting file

restart_dir = "restart"
restart_file = 1

rc_fname = "pgrid4py.dat" ## coord file to read PREFIX
d_fname =  "pygrid_tot.dat"    ## file to write PREFIX

MIN_RHO, FORCE_GEQ, MIN_PRS = 0,1,2
MIN_RHO_SHOCKPRS            = 4
####################################################################################################
###
### MIN_PRS:   Sets background to some minimum pressure, specified as a fraction of the maximum
### MIN_RHO:   Sets background to some minimum rho, specified as a fraction of the maximum
### FORCE_GEQ: Sets the background to the density minimum and picks P to force radial equilibrium
###
####################################################################################################


INIT_METHOD = MIN_RHO

####################################################################################################
###
### Set paramers for simluation
### Parameters for initial data (load from init file later, probably just easier to specify here)
###
####################################################################################################

NGAM = 3.       # polytropic index
MBH  = 2.e5     # mass of central object
MDISK= 2.e3     # total mass of disk
RDISK= 5       # radius of disk in units that we use
DIST = 1.04    # distortion parameter
RHOCRIT = 1.e-8 # minimum density as a fraction of maximum density
RHOCRIT_VISC = 1.e-1 ## same but for viscosity limiting
PRSCRIT = 1.e-8
GAMMA = 1 + 1/NGAM

### OTHER STUFF
GRAV_CONST = 6.67e-2  ## gravitational constants in the units we use
SPECAM     = np.sqrt(GRAV_CONST*MBH*RDISK)
## BASICALLY EVERYTHING BEFORE THIS DOESN'T REALLY MATTER

####################################################################################################
###
### Sets up the grid
###
####################################################################################################



#######################################################################
## need to do some stuff for parallelization 
#####################################################################

f = open(rc_fname, "r")
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
####################################################################################################
###
### defines that function that gives a (the closest) grid coordinate give a spatial point 
###
####################################################################################################

if restart_file:
    print("Using restart file")
    sim = pt.Simulation(restart_dir)
    rgr  = np.outer(rcs, np.ones(NTHCS) )
    thgr = np.outer(np.ones(NRCS), thcs)
    rhogr, prsgr, vx1gr, vx2gr, vx3gr = sim.gen_interpolate(rcs, thcs, nghosts=3)
    cgrid[:,:,X1] = rgr
    cgrid[:,:,X2] = thgr
    cgrid[:,:,RHO] = rhogr
    cgrid[:,:,PRS] = prsgr
    cgrid[:,:,VX1] = 0*vx1gr
    cgrid[:,:,VX2] = vx2gr
    cgrid[:,:,VX3] = vx3gr
else:
    print("Making clean initial conditions")
    ####################################################################################################
    ###
    ### Does the mass integral
    ###
    ####################################################################################################

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

    ####################################################################################################
    ###
    ### Sets up the mass profile (doesn't do anything about the ambient medium)
    ###
    ####################################################################################################

    def rho0(r, th):
        const = (GRAV_CONST * MBH/(NGAM+1)/entA/RDISK)**NGAM
        cdstuff = (RDISK/r-0.5*(RDISK/(r*np.sin(th)))**2-0.5/DIST)**NGAM
        return const*cdstuff
    ## this isn't the pressure that we use, but it's interesting to see how different
    ## it is from what actually produces an equilibrium state
    def pres0(rho):
        return entA*rho**GAMMA

    ## grav potential
    def phi(r):
        return -GRAV_CONST*MBH/r

    for i, r in enumerate(rcs):
        for j, th in enumerate(thcs):
            rho = rho0(r, th)
            cgrid[i,j,X1] = r
            cgrid[i,j,X2] = th
            if rho < 0:
                continue
            cgrid[i,j,RHO] = rho
            cgrid[i,j,PRS] = pres0(rho)
            cgrid[i,j,VX3] = SPECAM/( r*np.sin(th) )

    rho_max = np.max(cgrid[:,:,RHO])
    rho_crit= rho_max*RHOCRIT
    p0 = GRAV_CONST*MBH*rho_crit/rcs[-1]
    ## set the critical density
    for i, r in enumerate(rcs):
        for j, th in enumerate(thcs):
            if cgrid[i,j,RHO] < rho_crit:
                cgrid[i,j,RHO] = rho_crit
                cgrid[i,j,PRS] = 0

    min_prs = np.min([prs for prs in cgrid[:,:,PRS].flatten() if prs > 0])
    print("minimum pressure is ", min_prs)

    ############
    ##
    ## Theta boundary
    ##
    ############

    rb_index = -np.ones( (NTHCS, 2) )*np.inf
    for j, th in enumerate(thcs):
        for i, r in enumerate(rcs):
            if cgrid[i,j,RHO] == rho_crit:
                if i == NRCS - 1 or i == 0:
                    continue
                if cgrid[i+1,j,RHO]>rho_crit:
                    rb_index[j][0]=i
                elif cgrid[i-1,j,RHO]>rho_crit:
                    rb_index[j][1]=i

    ################
    ##
    ## Set up grid
    ##
    ################

    #for i, r in enumerate(rcs):
    #    for j, th in enumerate(thcs):
    #        if cgrid[i,j,RHO] == rho_crit:
    #            cgrid[i,j,PRS] = min_prs 
    #        #    cgrid[i,j,PRS] = min_prs #GRAV_CONST*MBH*rho_crit/r- p0

    ### background pressure grid in hydrostatic equilibrium
    ambprs_grid = np.zeros(np.shape(cgrid)[:-1])
    for i, r in enumerate(rcs):
        for j, th in enumerate(thcs):
            ambprs_grid[i][j] = GRAV_CONST*MBH*rho_crit/r

    if INIT_METHOD == MIN_RHO:
        min_ambprs = np.min(ambprs_grid)
        ambprs_grid -= min_ambprs - min_prs

    ### set up ambient pressure
    for i, r in enumerate(rcs):
        for j, th in enumerate(thcs):
            if cgrid[i,j,RHO] == rho_crit:
                cgrid[i,j,PRS] = ambprs_grid[i,j]
    ## check that we have the correct mass
    num_mass = 0
    for i, r in enumerate(rcs):
        for j, th in enumerate(thcs):
            rho = cgrid[i,j,RHO]
            if rho < rho_crit:
                continue
            if i == NRCS -1:
                continue
            if j == NTHCS - 1:
                continue

            dr = cgrid[i+1,j,X1] - cgrid[i,j,X1]
            dth= cgrid[i,j+1,X2] - cgrid[i,j,X2]
            vol_element = r**2*np.sin(th)*dr*dth
            num_mass += vol_element*rho

    num_mass *= 2*np.pi *2 #do the polar integral
    print(f"total mass fraction in grid is {num_mass/MDISK}")

    ####################################################################################################
    ###
    ### Sets up the ambient medium
    ###
    ####################################################################################################

    ## For the MIN_PRS AND MIN_RHO cases, we want pressures that have gradients ~GM*rho/r

    ## Find boundaries at each value of theta

    if INIT_METHOD == FORCE_GEQ:

        ## evolve from right boundary
        ## evolve from left boundary
        rb_index = -np.ones( (NTHCS, 2) )*np.inf
        for j, th in enumerate(thcs):
            for i, r in enumerate(rcs):
                if cgrid[i,j,RHO] == rho_crit:
                    if i == NRCS - 1 or i == 0:
                        continue
                    if cgrid[i+1,j,RHO]>rho_crit:
                        rb_index[j][0]=i
                    elif cgrid[i-1,j,RHO]>rho_crit:
                        rb_index[j][1]=i
        for j, th in enumerate(thcs):
            if rb_index[j][0] < 0 or rb_index[j][1] < 0:
                ## this means that there is not; set inner bou
                continue


            ## compute from left boundary
            rcl = rcs[:lindex][::-1]
            for i, r in enumerate(rcl):
                ieff = lindex - i -1
                rplus = cgrid[ieff+1,j,X1]
                dr = rplus - r
                rho = cgrid[ieff, j, RHO]

                gradph = (phi(rplus)-phi(r))/dr
                centf  = cgrid[ieff+1,j,VX3]**2/rplus*np.sin(th)
                dp = rho*dr*(centf-gradph)

                cgrid[ieff,j,PRS] = cgrid[ieff+1,j,PRS]-dp

            ## compute from right boundary
            rcr = rcs[rindex+1:]
            for i, r in enumerate(rcr):
                ieff = rindex+1+i
                if ieff == NRCS-1:
                    continue
                rminus = cgrid[ieff-1,j,X1]
                dr =  r - rminus

                rho = cgrid[ieff-1,j,RHO]

                gradph = (phi(r)-phi(rminus))/dr
                centf  = cgrid[ieff-1,j,VX3]**2/rplus*np.sin(th)
                dp = rho*dr*(centf-gradph)

                cgrid[ieff,j,PRS] = cgrid[ieff-1,j,PRS]+dp

        plboundary = GRAV_CONST*MBH*rho_crit/rcs[0]
        cgrid[:,:,PRS] += plboundary

        for j, th in enumerate(thcs):
            if rb_index[j][0] < 0 or rb_index[j][1] < 0:
                for i, r in enumerate(rcs):
                    rho = cgrid[i,j,RHO]
                    cgrid[i,j,PRS] = GRAV_CONST*MBH*rho/r

        cgrid[-1,:,PRS] = cgrid[-2,:,PRS]

    ## check if there are any negative pressures or densities
    print("minimum density", np.min(cgrid[:,:,RHO]))
    print("minimum pressure", np.min(cgrid[:,:,PRS]))

##the pressure gradients around the boundary are sometimes very big. 

####################################################################################################
###
### Writes to file
###
####################################################################################################


#########################################
###
### need to make pinit files for each process
###
###########################################


## input PCOORD for process -> output cgrid for pcoords
def write_cgrid():
    f   = open(d_fname, "w")
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
    print(f"Wrote file {d_fname}")

####################################################################################################
##
## write files
##
####################################################################################################

write_cgrid()

####################################################################################################
##
## Need to change the pluto.ini rhocrit parameter for viscous limiting
##
####################################################################################################

if False:
    ini_filer = open("pluto_template.ini", "r")
    lines = ini_filer.readlines()
    ini_filer.close()

    rhocrit_visc = RHOCRIT_VISC * rho_max

    for i,line in enumerate(lines):
        if line.split(" ")[0] == "RHOCRIT":
            print(f"WRITING RHOCRIT_VISC: {rhocrit_visc:5.5e}")
            lines[i] = f"RHOCRIT\t\t\t{rhocrit_visc:5.5e}"

            ini_filew = open("pluto.ini", "w")
            ini_filew.writelines(lines)
            ini_filew.close()
