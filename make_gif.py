import plutoplot as pp
import numpy as np
import matplotlib.pyplot as plt
import glob
from PIL import Image

sim = pp.Simulation("")
ts = sim.t
t_fin = 0
tindex = np.argmin(np.abs(ts - t_fin))
if t_fin==0:
    tindex = len(ts)-1
name = "test"

vx1_bd = np.zeros(2)
rho_bd = np.zeros(2)
## Get minimum and maximum pressures
for i in range(tindex):

    rho = sim[i].rho[:,-1,0]
    vx1 = sim[i].vx1[:,-1,0]
    
    if i == 0:
        vx1_bd = [np.min(vx1), np.max(vx1)]
        rho_bd = [np.min(rho), np.max(rho)]
        continue
    rho = np.append(rho, [rho_bd])
    vx1 = np.append(vx1, [vx1_bd])
    vx1_bd = [np.min(vx1), np.max(vx1)]
    rho_bd = [np.min(rho), np.max(rho)]    

## Make plots
for i in range(tindex):
    fig, ax1 = plt.subplots(1,1, figsize=(7,7))
    fig.suptitle(f"time: {ts[i]:03f}")
    ax1.set_xlabel(r"$r$")
    ax1.set_ylabel(r"$\rho$")
    
    ax2=ax1.twinx()
    
    ax2.set_ylabel(r"$v_r$")
    
    ax1.semilogy(sim[0].x1, sim[0].rho[:,-1,0], color="blue", linewidth=1, label=r"$v_{r,0}$")
    ax1.set_ylim(*rho_bd)
    ax1.semilogy(sim[-1].x1, sim[i].rho[:,-1,0], color="blue", linewidth=2, linestyle="--", label=r"$v_r$")

    ax2.set_xscale("linear")
    ax2.plot(sim[0].x1, sim[0].vx1[:,-1,0], color = "red", linewidth=1, label=r"$\rho_0$")
    ax2.set_ylim(*vx1_bd)
    ax2.plot(sim[-1].x1, sim[i].vx1[:,-1,0], color = "red", linewidth=2, linestyle="--", label=r"$\rho$")
    ax1.legend(loc=0, frameon=False)
    ax2.legend(loc=0, frameon=False)
    fig.savefig(f"imgs/{name}_{i:03d}.png", bbox_inches="tight", dpi=300)
    plt.close()

images = sorted(glob.glob("imgs/*.png"))
def make_gif(images):
    frames = [Image.open(image) for image in images]
    frame_one = frames[0]
    frame_one.save(f"{name}.gif", format="GIF", append_images=frames,
               save_all=True, duration=100, loop=0)

make_gif(images)
