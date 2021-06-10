using PyPlot
using FFTW
using Plots

include("Poisson_utils.jl")
include("Poisson_const.jl")
PyPlot.close("all")
## Poisson Process

simtime = 300.0
times = collect(0.0:dt:simtime)
steps = length(times)

params = [sig./freq, freq, gain, offset]
raster, ggs = osc_poisson(dt, simtime, steps, params, g=true,dephase = 0.0)
raster2, ggs2 = osc_poisson(dt, simtime, steps, params, g=true,dephase = 0.25)
raster3, ggs3 = osc_poisson(dt, simtime, steps, params, g=true,dephase = 0.5)
raster4, ggs4 = osc_poisson(dt, simtime, steps, params, g=true,dephase = 0.75)

my_dpi=96
fig2,ax2 = PyPlot.subplots(4,2,figsize=(900/my_dpi, 600/my_dpi), dpi=my_dpi)
ax2[1,1].plot(times,ggs)
ax2[1,1].axis("off")
ax2[1,1].set_xlim(100,300)
ax2[1,2].plot(times,raster,c="k")
ax2[1,2].axis("off")
ax2[1,2].set_xlim(100,300)

ax2[2,1].plot(times,ggs2,c="orange")
ax2[2,1].set_xlim(100,300)
ax2[2,1].axis("off")
ax2[2,2].plot(times,raster2,c="k")
ax2[2,2].axis("off")
ax2[2,2].set_xlim(100,300)

ax2[3,1].plot(times,ggs3,c="green")
ax2[3,1].axis("off")
ax2[3,1].set_xlim(100,300)
ax2[3,2].plot(times,raster3,c="k")
ax2[3,2].axis("off")
ax2[3,2].set_xlim(100,300)

ax2[4,1].plot(times,ggs4,c="red")
ax2[4,1].axis("off")
ax2[4,1].set_xlim(100,300)
ax2[4,2].plot(times,raster4,c="k")
ax2[4,2].axis("off")
ax2[4,2].set_xlim(100,300)
