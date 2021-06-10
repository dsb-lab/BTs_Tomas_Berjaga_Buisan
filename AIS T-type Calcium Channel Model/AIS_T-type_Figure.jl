#Import all the packages
using Pkg
using PyCall
using PyPlot
pyplt=PyPlot;
my_dpi= 96;

include("network_utils.jl")
include("AIS_T-type_const.jl")
include("AIS_T-type_utils.jl")

#Activation and Inactivation Gating Variables
voltage=collect(-120.0:0.1:0.0);
m,h = gating_variables(m_med_t,h_med_t,k1,k2,voltage)
fig, ax = pyplt.subplots(figsize=(700/my_dpi, 500/my_dpi), dpi=my_dpi)
ax.plot(voltage,m,"k")
ax.set_xlabel("Voltage (mV)", fontsize=20)
ax.tick_params(axis="x", labelsize=15)
ax.set_ylabel("Activation", fontsize=20)
ax.set_xlim(-130,10)
ax.tick_params(axis="y", labelsize=10)
ax.spines["top"].set_visible(false) # Hide the top edge of the axis
ax.spines["right"].set_visible(false) # Hide the right edge of the axis
ax.spines["left"].set_position(("axes",-0.02)) # Offset the left scale from the axis
ax.spines["bottom"].set_position(("axes",-0.02)) # Offset the bottom scale from the axis

ax2 = ax.twinx()
ax2.plot(voltage,h,"k")
ax2.set_xlabel("Voltage (mV)", fontsize=20)
ax2.tick_params(axis="x", labelsize=15)
ax2.set_ylabel("Inactivation", fontsize=20)
ax2.tick_params(axis="y", labelsize=10)
ax2.spines["top"].set_visible(false) # Hide the top edge of the axis
ax2.spines["left"].set_position(("axes",-0.02)) # Offset the left scale from the axis
ax2.spines["bottom"].set_position(("axes",-0.02)) # Offset the bottom scale from the axis
fig.tight_layout()

# tau_m and tau_h
voltage=collect(-120:3.0:0);
tau_m_1 = compute_tau_m(voltage)
fig1, ax1 = pyplt.subplots(figsize=(900/my_dpi, 600/my_dpi), dpi=my_dpi)
ax1.plot(voltage,tau_m_1,"k",voltage,tau_m_1,"ks")
ax1.set_xlabel("Voltage (mV)", fontsize=20)
ax1.tick_params(axis="x", labelsize=15)
ax1.set_ylabel("\$\\tau_m\$ (ms)", fontsize=20)
ax1.set_xlim(-130,10)
ax1.tick_params(axis="y", labelsize=15)
ax1.spines["top"].set_visible(false) # Hide the top edge of the axis
ax1.spines["right"].set_visible(false) # Hide the right edge of the axis
ax1.spines["left"].set_position(("axes",-0.02)) # Offset the left scale from the axis
ax1.spines["bottom"].set_position(("axes",-0.02)) # Offset the bottom scale from the axis
fig1.tight_layout()

voltage1=collect(-120:0.01:-75);
voltage2=collect(-70:0.01:0);
tau_h_1 = compute_tau_h(voltage1)
tau_h_2 = compute_tau_h(voltage2)

fig2, ax2 = pyplt.subplots(figsize=(900/my_dpi, 600/my_dpi), dpi=my_dpi)
ax2.plot(voltage1,tau_h_1,c="k",label=" \$\\tau_{rec}\$")
ax2.set_ylim(0,500)
ax2.set_xlabel("Voltage (mV)", fontsize=20)
ax2.set_ylabel(" \$\\tau_{rec}\$ (ms)", fontsize=20)
ax2.tick_params(axis="y", labelsize=10)
ax2.tick_params(axis="x", labelsize=15)
ax2.spines["top"].set_visible(false) # Hide the top edge of the axis
ax2.spines["left"].set_position(("axes",-0.02)) # Offset the left scale from the axis
ax2.spines["bottom"].set_position(("axes",-0.02)) # Offset the bottom scale from the axis

ax3 = ax2.twinx()
ax3.plot(voltage2,tau_h_2,"b",label="\$\\tau_h\$")
ax3.set_xlabel("Voltage (mV)", fontsize=20)
ax3.tick_params(axis="x", labelsize=15)
ax3.set_ylabel("\$\\tau_h\$ (ms)", fontsize=20)
ax3.set_ylim(0,150)
ax3.tick_params(axis="y", labelsize=10)
ax3.spines["top"].set_visible(false) # Hide the top edge of the axis
ax3.spines["left"].set_position(("axes",-0.02)) # Offset the left scale from the axis
ax3.spines["bottom"].set_position(("axes",-0.02)) # Offset the bottom scale from the axis

fig2.legend(loc="upper right", ncol=3, borderaxespad=4,fontsize=20)
fig2.tight_layout()

#Maximal Current per stationary membrane potential
V = collect(-100:10:50)
tvec = collect(0:0.1:2000)
I = zeros(length(V))
for (i,t) in enumerate(V)
    m_t= 1.0 / (1.0+exp(-(t-m_med_t)/k1))
    h_t= 1.0 / (1.0+exp((t-h_med_t)/k2))

    I[i] = 75.0*(m_t^2)*h_t*(120-t)
end

fig4,ax4 = pyplt.subplots(figsize=(700/my_dpi, 500/my_dpi), dpi=my_dpi)
ax4.plot(V,I,c="k")
ax4.invert_yaxis()
ax4.set_ylabel("Current (pA)",loc="bottom", fontsize=15)
ax4.set_xlabel("Voltage (mV)", loc="left", fontsize=15)
ax4.spines["top"].set_visible(false) # Hide the top edge of the axis
ax4.spines["right"].set_visible(false)
ax4.spines["bottom"].set_position("zero")
ax4.spines["left"].set_position("center")
ax.tick_params(axis="x", labelsize=15)
ax.tick_params(axis="y", labelsize=15)
fig4.tight_layout()
