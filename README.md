# bnsm-disk
some kind of documentation

## definitions.h and pluto.ini
### Length scales and variables
In order to add in approximate r-process heating rates, you need to specify physical units. The three relevant ones are

`UNIT_DENSITY` $\rho = 10^8 \text{g }\text{cm}^{-3}$

`UNIT_LENGTH` $L = 10^7 \text{cm }$

`UNIT_VELOCITY` $v = 10^8 \text{cm }\text{s}^{-1}$

Next unit scales are not fundamental, but are good to calculate/know

`UNIT_MASS` $M = L^3\rho = 10^{29} \text{g} \simeq 5\times 10^{-5} M_\odot$

`UNIT_TIME` $t = L/v = 0.1 \text{s}$

`UNIT_ENERGY` $U = Mv^2=10^{45} \text{ergs}$

Really the only fundamental constant that we need to convert is $G$. Note that $G\sim v^2L/M$ so we should scale by the inverse of this, i.e.

`CONST_G` $G = 6.67\times 10^{-8} \text{ (cgs units)} \times \simeq 6.67\times 10^{-2} \text{ (code units)}$

Probably the most important quantity is the r-process heating rate $q\sim E/(Mt)$, meaning that we should take for our code-unit heating rate

`RP_HEATING` $q = 4\times 10^{16}\text{ }t^{-1.3} \text{cgs units} = 0.4\text{ }t^{-1.3}\text{ (code units)}.$

### User defined parameters

There are currently 7 user defined parameters, and I think all are currently used.

`MBH` Mass of central object. We use $M_{BH} = 2\times 10^5 \simeq 0.2\text{ }M_\odot$

`NGAM` Polytopic index. We take $n=3$ corresponding to $\gamma=4/3$ for a radiation-dominated gas,

`TSTART` Time (in code units) when the viscosity is turned on. 

`ALPHA` Usual $\alpha$ viscosity parameter

`VISCLIM` Whether or not to limit the viscosity is set, according to whatever procedure is set (now done according to a combination of tracers/ad-hoc source limiting which I want to remove). 

`RHOCRIT` Sets the critical density of the disk where the tracers are initially set.

`EXP_LIM` Sets how rapidly the viscosity source terms are set to fall off outside of the region tracer particles have access to.

The mass of the disk is probably set in `initgrid.py` and disn information. 

## init.c

`Init` Initializes the tracers based on `RHOCRIT.` Initialization of the density and velocity data is done in `main.c` and `startup.c`

`BodyForcePotential` Set central potential $\Phi=-G M_\text{BH}/r$

`UserDefBoundary` Had some issues with mass inflow at the outer radial boundary and inner polar boundary. Imposed no inflow boundary conditions there, which are almost identical to the usual Neumann ones except we set $v_r$ (or $v_\phi$) to zero if they could represent inflow. 

## initgrid.py and initgridn.py

These are used to initialize the grid. I don't know enough C well

## main.c



## radiat.c

## startup.c

## visc_nu.c

## viscous_rhs.c
