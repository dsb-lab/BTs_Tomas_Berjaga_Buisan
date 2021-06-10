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

##
# This code is prepared for simulating and saving the Network Simulations for
# future analisis. The spike step, index, a random pyramidal trace, matrix connectivity and Post Synaptic Conductance
# The only change that is needed to be functional is to change the Directories.
##
include("network_constants.jl")
include("network_utils.jl")
include("PN_SOM_PV_solver_t_type.jl")

#Time parameters for the simulations
    Tmax=3000.0 ;
    dt = 0.1 ;
    const delay_time=1.0 ;
    const Dmax = delay_time/dt;
    tvec = collect(0.0:dt:Tmax);

#Number of entities in the network
    number_pyramidal=800;
    number_PV=100;
    number_SOM=100;
    number_total=number_pyramidal+number_PV+number_SOM;

#Connectivity Matrix Probabilities
    #PYR NETWORK (HERE YOU CAN CHANGE FOR ANY OF THE POSSIBLE MATRICES)
    prob_connect_pyr_pyr=0.1 ; prob_connect_pyr_PV=0.00;  prob_connect_pyr_SOM=0.0;
    prob_connect_PV_PV=0.00; prob_connect_PV_pyr=0.0; prob_connect_PV_SOM=0.0;
    prob_connect_SOM_SOM=0.00 ; prob_connect_SOM_pyr=0.00 ; prob_connect_SOM_PV=0.0;
    #PV
    I_inject_PV_initial = 100.0 ;
    #SOM
    I_inject_SOM_initial = 100.0 ;
    #AIS T-type Parameters
    const E_Ca_t = 120.0;  #Mirar
    const k1 = 5.2;
    const k2 = 6.2;
    const m_med_t = -48.4; #CaV 3.2
    const h_med_t = -75.6;

#Parameters packing
par_number = (number_pyramidal,number_PV,number_SOM);
par_connectivity = (prob_connect_pyr_pyr,prob_connect_pyr_PV,prob_connect_pyr_SOM,prob_connect_PV_PV,prob_connect_PV_pyr, prob_connect_PV_SOM,prob_connect_SOM_SOM,prob_connect_SOM_pyr,prob_connect_SOM_PV);

par_time=(tvec,dt,Dmax,delay_time);

#=
g_Ca_t = 0.0;

I_inject_S_initial=750.0;I_inject_D_initial=750.0;
par_initial=(I_inject_S_initial,I_inject_D_initial,I_inject_PV_initial,I_inject_SOM_initial);

#Solving Net
for simulations in 1:10
    Matrix_connectivity = connectivity_network(par_number,par_connectivity);
    @time steps_pyr,index_pyr,steps_all,index_all,V_S,V_D= PN_SOM_PV_solve_T_type(par_time,par_number,par_initial,g_Ca_t,Matrix_connectivity);
    S = PSC(dt, steps_pyr,tvec);
    S = S[5000:end]
    tvec2 = tvec[5000:end]
    writedlm("PYR NETWORK/g=0/I=750/PSC_t_Inject_750_PYR_$simulations.txt", S)
    writedlm("PYR NETWORK/g=0/I=750/connectivity/Matrix_Inject_750_PYR_$simulations.txt", Matrix_connectivity)
    writedlm("PYR NETWORK/g=0/I=750/steps/Spike_steps_Inject_750_PYR_$simulations.txt", steps_all)
    writedlm("PYR NETWORK/g=0/I=750/index/Spike_index_Inject_750_PYR_$simulations.txt", index_all)
    writedlm("PYR NETWORK/g=0/I=750/voltages/VS_Inject_750_PYR_$simulations.txt", V_S)
    writedlm("PYR NETWORK/g=0/I=750/voltages/VD_Inject_750_PYR_$simulations.txt", V_D)
end
=#
