# bnsm-disk
some kind of documentation

## definitions.h
### Length scales and variables
In order to add in approximate r-process heating rates, you need to specify physical units. The three relevant ones are

`UNIT_DENSITY` $\rho = 10^8 \text{g }\text{cm}^{-3}$

`UNIT_LENGTH` $L = 10^7 \text{cm }$

`UNIT_VELOCITY` $v = 10^8 \text{cm }\text{s}^{-1}$

Next unit scales are not fundamental, but are good to calculate/know

`UNIT_MASS` $M = L^3\rho = 10^{29} \text{g} \simeq 5\times 10^{5} M_\odot$

`UNIT_TIME` $t = L/v = 0.1 \text{s}$

`UNIT_ENERGY` $U = Mv^2=10^{45} \text{ergs}$

Really the only fundamental constant that we need to convert is $G$. Note that 

$$G\sim v^2L/M$$

so we should scale by the inverse of this, i.e.

`CONST_G` $G = 6.67\times 10^{-8} \text{ (cgs units)} \times \simeq 6.67\times 10^{-2} \text{ (code units)}$

Probably the most important quantity is the r-process heating rate $q\sim E/(Mt)$, meaning that we should take for our code-unit heating rate

`RP_HEATING` $q = 4\times 10^{16}\text{ }t^{-1.3} \text{cgs units} = 0.4\text{ }t^{-1.3}\text{ (code units)}.$

### User defined parameters

## init.c

## initgrid.py and initgridn.py

## main.c

## pluto.ini

## radiat.c

## startup.c

## visc_nu.c

## viscous_rhs.c
