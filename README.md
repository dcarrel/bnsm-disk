# bnsm-disk
some kind of documentation

## definitions.h and pluto.ini
### Length scales and variables
In order to add in approximate r-process heating rates, one needs to specify physical units. The three relevant ones are

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

Note that in the above $t$ has to be divided by the appropriate `UNIT_TIME.`

### User defined parameters

There are currently 7 user defined parameters, and I think all are currently used.

`MBH` Mass of central object. We use $M_{BH} = 2\times 10^5 \simeq 0.2\text{ }M_\odot$

`NGAM` Polytopic index. We take $n=3$ corresponding to $\gamma=4/3$ for a radiation-dominated gas,

`TSTART` Time (in code units) when the viscosity is turned on. 

`ALPHA` Usual $\alpha$ viscosity parameter

`VISCLIM` Whether or not to limit the viscosity is set, according to whatever procedure is set (now done according to a combination of tracers/ad-hoc source limiting which I want to remove). Don't think this is actually used anymore.

`RHOCRIT` Sets the critical density of the disk where the tracers are initially set.

`EXP_LIM` Sets how _small_ (used to be how large) the ratio between the conservative variable and source has to be in order to do the explicit limiting.

The mass of the disk is probably set in `initgrid.py` and disn information. 

## init.c

`Init` Initializes the tracers based on `RHOCRIT.` Initialization of the density and velocity data is done in `main.c` and `startup.c`

`BodyForcePotential` Set central potential $\Phi=-G M_\text{BH}/r$

`UserDefBoundary` Had some issues with mass inflow at the outer radial boundary and inner polar boundary. Imposed no inflow boundary conditions there, which are almost identical to the usual Neumann ones except we set $v_r$ (or $v_\phi$) to zero if they could represent inflow. 

## set_grid.c/startup.c/ and initgrid.py/initgridn.py

It's easier (for me) to set up and load the initial data with Python. The way this works is as follows:

1. `set_grid.c` defines the global grid, and writes the grid to `pgrid4py.dat.'
2. A system call is then made to `initgrid.py` which will write the initial primitives to `pygrid_tot.dat.` This is probably the most important step. See [here](https://github.com/dcarrel/bnsm-disk/blob/main/acc_disk_hydro.pdf) for the details of disk initialization. 
3. MPI assigns different grid regions to different processesors and so to initialize the generated data, one has to do this subgrid by subgrid. This is done in `startup.c.` Each process does a system call to `initgridn.py.` This generates process-specific grid files `pygrid_nnn.dat` which is written to the initial PLUTO grid.

## radiat.c

This is where the heating rate is implemented. To use a given functional form, one needs to set the cooling rate to `TABULATED` in `definitions.h` and then define `rhs[RHOE]` appropriately (positive for heating). 

## visc_nu.c

Here, we use Sunyaev-Shakura $\alpha$ viscosity prescription: 

$$\nu = \frac{\alpha P}{\rho r^{3/2}\sqrt{G M_{BH}}}$$

This is set in `visc_nu.c` which requires that `VISCOSITY` be set to `EXPLICIT` in `definitions.h`. In regions where the sounds speed is high (like the ambient medium around the disk), the viscous source terms are very large and cause numerical problems. There are two ways to deal with this

1. Limit the viscous source terms only to regions where the disk is spread. This can be done by using tracers, which are an _almost_ binary parameter that is `1` in places where the disk has advected to and `0` elsewhere.
2. Explicitly limit viscosity by reducing the source term where it is too large. This is done in `viscous_rhs.c,` and is used in tandem with (1). 

## viscous_rhs.c

The viscous source is limited via `LIMIT_RHS,` which has not a great selection criteria right now. I will change it and update this.

There is also a user-defined parameter `viscon.` This is set concurrently with (but not by) calls to `LIMIT_RHS` and is `1` where the viscous source is one and `0` where it is off. 
