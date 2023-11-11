import plutoplot as pp
import numpy as np
import matplotlib.pyplot as plt

sim = pp.Simulation("")
siml = len(sim)
x1s, x2s = sim[0].x1, sim[0].x2

x2s_index = len(x2s) - np.logspace(0, np.log10(len(x2s)), 10)
x2s_index = np.array(np.floor(x2s_index)[:-1], dtype="int")

RHO, PRS, VX1, VX2, VX3 = range(0,5)
PNAMES = ["rho", "prs", "vx1", "vx2", "vx3"]
CS, MN = range(0,2)

min_vals, max_vals = np.zeros(6), np.zeros(6)

def update_pext(prim, thslice):
    mi,ma = np.min(thslice), np.max(thslice)
    if min_vals[prim] > mi:
        min_vals[prim] = mi
    if max_vals[prim] < ma:
        max_vals[prim] = ma
pgrid = np.zeros((len(x2s_index),5,len(x1s)))
vgrid = np.zeros((len(x2s_index),2,len(x1s)))
    
for j, s in enumerate(sim):
  for i, x2i in enumerate(x2s_index):
      th = x2s[x2i]
      ##update prim extrema
      for p in range(0,5):
          

      for k, r in enumerate(x1s):  
          if True:
              vx1 = s[PNAMES[VX1]][k, x2i,0]
              vx2 = s[PNAMES[VX2]][k,x2i,0]
              vx3 = s[PNAMES[VX3]][k,x2i,0]
              vmag = np.sqrt(vx1**2 + r**2*vx2**2
                            + r**2*np.sin(th)**2*vx3**2)

              rho = s[PNAMES[RHO]][k, x2i, 0]
              prs = s[PNAMES[PRS]][k, x2i, 0]
              #rhor = s[PNAMES[RHO]][k:k+2, x2i,0] 
              #prsr = s[PNAMES[PRS]][k:k+2, x2i,0]
              #rhoth = s[PNAMES[RHO]][k,x2i:x2i-2,0]
              #prsth = s[PNAMES[PRS]][k,x2i:x2i-2,0]
              vinf[i][0][k]= np.sqrt(prs/rho)#np.sqrt((prsr[1]-prsr[0])/(rhor[1]-rhor[0])
                         #+ (prsth[1]-prsth[0]) 
              vinf[i][1][k] = vmag/np.sqrt(prs/rho)
      
