import plutoplot as pp
import numpy as np
import sys
from datetime import datetime
import matplotlib.pyplot as plt

nsample = 5 ## sample every five
tstart  = datetime.today().strftime("%Y%m%d_%H%M%S")
print(tstart)
try:
    nsample = sys.argv[1]
except:
    print("using default sample")
    
## remove old images
    
sim = pp.Simulation("")
nsteps = len(sim)
sidxs = np.arange(0, nsteps, nsample)

rcs  = sim.x1
thcs = sim.x2 

prims = range(0,5)
RHO,PRS,VX1,VX2,VX3 = prims
names = ["rho","prs","vx1","vx2","vx3"]

plot_grid = False

## make density plots
if plot_grid:
    for i in sidxs:
        simt = sim[i]
        for p in prims:
            ft, at = simt.plot(names[p])
            ft.savefig(f"imgs/{names[p]}/{tstart}_{names[p]}_{i:03d}.png", bbox_inches="tight")
            plt.close(ft)

## Better plots are probably just radial lines
rhocrit = 1e-2
numlines= 10
for i in sidxs:
    simt = sim[i]
    rts  = simt.x1
    thts = simt.x2
    parr = np.zeros((5,len(rts), len(thts)))
    for p in prims:
        parr[p] = simt[names[p]][:,:,0]

    rhomin = rhocrit*np.max(parr[RHO])
    dtharray = (parr[RHO] - rhomin).T
    thmin = np.min( [i for i, dth in enumerate(dtharray) if np.max(dth) > 0] )
    thindx_samp = np.array(len(thts)-np.logspace(0, np.log10(len(thts)), numlines), dtype="int")
    ## find minimum and maximum values
    min_vals = np.zeros(5)
    max_vals = np.zeros(5)
    for thindex in thindx_samp:
        for p in prims:
            cmax = np.max(parr[p].T[thindex])
            cmin = np.min(parr[p].T[thindex])
            if cmax > max_vals[p]:
                max_vals[p] = cmax
            if cmin < min_vals[p]:
                min_vals[p] = cmin
                
    ## actually do the plotting
    for p in prims:
      fig, axs = plt.subplots(1,1, figsize=(5,5))
      mint, maxt = min_vals[p], max_vals[p]
      if not (mint < maxt):
          mint, maxt = -1,1
      axs.set_ylim(mint, maxt)
      axs.set_xlabel("r")
      axs.set_ylabel("th")
      cmap = plt.cm.get_cmap("viridis")
      for j,thindex in enumerate(thindx_samp):
          axs.semilogx(rts, parr[p].T[thindex], color = cmap(j/len(thindx_samp)))
      fig.suptitle(f"t={simt.t}, {names[p]}")
      fig.savefig(f"imgs/{names[p]}/{tstart}_lp{names[p]}_{i:03d}.png", bbox_inches="tight")
      plt.close(fig)
            

    


    
