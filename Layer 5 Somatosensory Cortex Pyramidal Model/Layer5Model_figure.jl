#Import all the packages
using Pkg
using Plots
using Parameters
using Random, Distributions
using FFTW
using PyCall
using PyPlot
using Profile
using DelimitedFiles
pyplt=PyPlot;

include("Poisson_utils.jl")
include("Poisson_const.jl")
include("Layer5Model_const.jl")
include("Layer5Model_utils.jl")
include("Spike_Bursts_utils.jl")

Tmax=3000.0 ;
dt = 0.1 ;
const delay_time=1.0 ;
const D=0 ;
const Dmax = delay_time/dt;
tvec = collect(0.0:dt:Tmax);
t_length=length(tvec)

freq = 10
sig= 100.0
gain= 30.0
offset = 0.0
params = [sig/freq, freq, gain, offset]
ts=[]
for i=1:50
    raster, ggs = osc_poisson(dt, Tmax, t_length , params, g=true,dephase = rand([0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4]))
    append!(ts, compute_ts(raster))
end
poisson_I= PSC(dt, ts, tvec)
Inject_D=0

soma,dendrite,T,poisson_current,II_Ca,II_AHP=pyr_solver(tvec,dt,poisson_I,Inject_D,0)
b_spks_150_9, b_ev_150_9 = burst_detection(soma,5.0)

my_dpi=96
fig2,ax2 = PyPlot.subplots(4,1,figsize=(1200/my_dpi, 800/my_dpi), dpi=my_dpi)
ax2[1,1].plot(tvec,soma,c="k",label="\$V_S\$")
ax2[1,1].plot(tvec,dendrite,c="c",label="\$V_D\$")
ax2[1,1].set_xlim(1000,2000)
ax2[1,1].set_ylabel("Voltage (mV)",fontsize=15)
ax2[1,1].legend(loc="best")
ax2[1,1].spines["top"].set_visible(false) # Hide the top edge of the axis
ax2[1,1].spines["right"].set_visible(false)

ax2[3,1].plot(tvec,II_AHP,c="k",label="\$I_{AHP}\$")
ax2[3,1].plot(tvec,II_Ca,c="c",label="\$I_{HVA}\$")
ax2[3,1].set_xlim(1000,2000)
ax2[3,1].legend(loc="upper right")
ax2[3,1].set_ylabel("Current (pA)",fontsize=15)
ax2[3,1].spines["top"].set_visible(false) # Hide the top edge of the axis
ax2[3,1].spines["right"].set_visible(false)

ax2[4,1].plot(tvec,T,c="k",label="AIS T-type current")
ax2[4,1].set_xlim(1000,2000)
ax2[4,1].set_ylim(-0.05,0.05)
ax2[4,1].legend(loc="best")
ax2[4,1].set_ylabel("Current (pA)",fontsize=15)
ax2[4,1].set_xlabel("Time (ms)",fontsize=15)
ax2[4,1].spines["top"].set_visible(false) # Hide the top edge of the axis
ax2[4,1].spines["right"].set_visible(false)

ax2[2,1].plot(tvec,poisson_current,c="k",label="Poisson current")
ax2[2,1].set_xlim(1000,2000)
ax2[2,1].legend(loc="upper right")
ax2[2,1].set_ylabel("Current (pA)",fontsize=15)
ax2[2,1].spines["top"].set_visible(false) # Hide the top edge of the axis
ax2[2,1].spines["right"].set_visible(false)
fig2.tight_layout()
