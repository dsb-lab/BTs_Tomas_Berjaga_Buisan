#Import all the packages
using Pkg
using PyCall
using PyPlot
using Distributions
pyplt=PyPlot;
include("PV_utils.jl")
include("PV_const.jl")
include("network_utils.jl")

my_dpi= 96;
Tmax, dt  = 2000, 0.01 ; #ms
tvec = collect(0:dt:Tmax);

v = zeros(length(tvec));
u = zeros(length(tvec));
v[1]=-55;

I=step_current(dt,tvec,50,150,300)
v,u,I,spikes = simulate_PV(tvec,v,u,I)

fig, ax = pyplt.subplots(2,1,figsize=(1000/my_dpi, 900/my_dpi), dpi=my_dpi)
ax[1].plot(tvec,v,"black")
ax[1].tick_params(axis="x", labelsize=15)
ax[1].tick_params(axis="y", labelsize=10)
ax[1].set_ylabel("Voltage (mV)", fontsize=20)
ax[1].set_xlim(0,200)
ax[1].spines["top"].set_visible(false) # Hide the top edge of the axis
ax[1].spines["right"].set_visible(false) # Hide the right edge of the axis

ax[2].plot(tvec,I,"black")
ax[2].set_xlabel("Time (ms)", fontsize=20)
ax[2].set_xlim(0,200)
ax[2].tick_params(axis="x", labelsize=15)
ax[2].set_ylim(0,400)
ax[2].tick_params(axis="y", labelsize=10)
ax[2].set_ylabel("Current (pA)", fontsize=20)
ax[2].spines["top"].set_visible(false) # Hide the top edge of the axis
ax[2].spines["right"].set_visible(false) # Hide the right edge of the axis
fig.tight_layout()

#Frequency Plot
v = zeros(length(tvec));
u = zeros(length(tvec));
v[1]=-55;
I = collect(0:25:500)
@time line_frequencies=frequency(dt,tvec,v,u,I)
frequencias=convert(Vector{Int64}, line_frequencies)
fig1, ax1 = pyplt.subplots(figsize=(1200/my_dpi, 600/my_dpi), dpi=my_dpi)
ax1.plot(I,frequencias,"black")
ax1.spines["top"].set_visible(false) # Hide the top edge of the axis
ax1.spines["right"].set_visible(false) # Hide the right edge of the axis
ax1.set_xlabel("Injected Current (pA)", fontsize=20)
ax1.tick_params(axis="x", labelsize=15)
ax1.tick_params(axis="y", labelsize=10)
ax1.set_ylabel("Frequency (Hz)", fontsize=20)

#Ornstein-Uhlenbeck process
v = zeros(length(tvec));
u = zeros(length(tvec));
v[1]=-55;

const sigma=300.0; const tau =3.0;
const i_tau_2 = 2.0/tau;
I_Inject_PV=zeros(length(tvec))
I_Inject_PV[1]=100;
for (i,t) in enumerate(tvec[1:end-1])
    I_Inject_PV[i+1] = I_Inject_PV[i] + ((100-I_Inject_PV[i])/ tau)*dt + sigma*rand(Normal(0, 1))*sqrt(i_tau_2*dt);
end

v,u,I,spikes = simulate_PV(tvec,v,u,I_Inject_PV)
fig2, ax2 = pyplt.subplots(2,1,figsize=(1000/my_dpi, 900/my_dpi), dpi=my_dpi)
ax2[1].plot(tvec,v,"black")
#ax[1].set_xlabel("Time (ms)", fontsize=15)
ax2[1].tick_params(axis="x", labelsize=15)
ax2[1].tick_params(axis="y", labelsize=15)
ax2[1].set_ylabel("Voltage (mV)", fontsize=20)
ax2[1].set_xlim(0,200)
ax2[1].spines["top"].set_visible(false) # Hide the top edge of the axis
ax2[1].spines["right"].set_visible(false) # Hide the right edge of the axis


ax2[2].plot(tvec,I,"black")
ax2[2].set_xlabel("Time (ms)", fontsize=20)
ax2[2].tick_params(axis="x", labelsize=15)
ax2[2].tick_params(axis="y", labelsize=15)
ax2[2].set_ylabel("Current (pA)", fontsize=20)
ax2[2].set_xlim(0,200)
ax2[2].spines["top"].set_visible(false) # Hide the top edge of the axis
ax2[2].spines["right"].set_visible(false) # Hide the right edge of the axis
fig2.tight_layout()
