using Pkg
using PyCall
using PyPlot
using Profile
using DelimitedFiles
pyplt=PyPlot;

include("Spike_Bursts_utils.jl")
##

#Please consider that to run this file, you must have first saved some network simulations.
#The file is prepared for the analysis of burst events for 10 networks for different AIS T-type Conductances.
#The only change that is needed to be functional is to change the Directories.

##
#=
##g=0
steps1 = readdlm("SOM NETWORK/g=0/I=250/steps/Spike_steps_Inject_250_SOM_1.txt",Int64)
index1 = readdlm("SOM NETWORK/g=0/I=250/index/Spike_index_Inject_250_SOM_1.txt",Int64)

steps2 = readdlm("SOM NETWORK/g=0/I=250/steps/Spike_steps_Inject_250_SOM_2.txt",Int64)
index2 = readdlm("SOM NETWORK/g=0/I=250/index/Spike_index_Inject_250_SOM_2.txt",Int64)

steps3 = readdlm("SOM NETWORK/g=0/I=250/steps/Spike_steps_Inject_250_SOM_3.txt",Int64)
index3 = readdlm("SOM NETWORK/g=0/I=250/index/Spike_index_Inject_250_SOM_3.txt",Int64)

steps4 = readdlm("SOM NETWORK/g=0/I=250/steps/Spike_steps_Inject_250_SOM_4.txt",Int64)
index4 = readdlm("SOM NETWORK/g=0/I=250/index/Spike_index_Inject_250_SOM_4.txt",Int64)

steps5 = readdlm("SOM NETWORK/g=0/I=250/steps/Spike_steps_Inject_250_SOM_5.txt",Int64)
index5 = readdlm("SOM NETWORK/g=0/I=250/index/Spike_index_Inject_250_SOM_5.txt",Int64)

steps6 = readdlm("SOM NETWORK/g=0/I=250/steps/Spike_steps_Inject_250_SOM_6.txt",Int64)
index6 = readdlm("SOM NETWORK/g=0/I=250/index/Spike_index_Inject_250_SOM_6.txt",Int64)

steps7 = readdlm("SOM NETWORK/g=0/I=250/steps/Spike_steps_Inject_250_SOM_7.txt",Int64)
index7 = readdlm("SOM NETWORK/g=0/I=250/index/Spike_index_Inject_250_SOM_7.txt",Int64)

steps8 = readdlm("SOM NETWORK/g=0/I=250/steps/Spike_steps_Inject_250_SOM_8.txt",Int64)
index8 = readdlm("SOM NETWORK/g=0/I=250/index/Spike_index_Inject_250_SOM_8.txt",Int64)

steps9 = readdlm("SOM NETWORK/g=0/I=250/steps/Spike_steps_Inject_250_SOM_9.txt",Int64)
index9 = readdlm("SOM NETWORK/g=0/I=250/index/Spike_index_Inject_250_SOM_9.txt",Int64)

steps10 = readdlm("SOM NETWORK/g=0/I=250/steps/Spike_steps_Inject_250_SOM_10.txt",Int64)
index10 = readdlm("SOM NETWORK/g=0/I=250/index/Spike_index_Inject_250_SOM_10.txt",Int64)

g_0_1,g_0_2,g_0_3,g_0_4,g_0_5,g_0_6,g_0_7,g_0_8,g_0_9,g_0_10=spike_events_all(steps1,index1,steps2,index2,steps3,index3,steps4,index4,steps5,index5,steps6,index6,steps7,index7,steps8,index8,steps9,index9,steps10,index10)

mean_b_spks_0_1=0;
mean_b_ev_0_1 =0;

mean_b_spks_0_2=0;
mean_b_ev_0_2 =0;

mean_b_spks_0_3=0;
mean_b_ev_0_3 =0;

mean_b_spks_0_4=0;
mean_b_ev_0_4 =0;

mean_b_spks_0_5=0;
mean_b_ev_0_5 =0;

mean_b_spks_0_6=0;
mean_b_ev_0_6 =0;

mean_b_spks_0_7=0;
mean_b_ev_0_7 =0;

mean_b_spks_0_8=0;
mean_b_ev_0_8 =0;

mean_b_spks_0_9=0;
mean_b_ev_0_9 =0;

mean_b_spks_0_10=0;
mean_b_ev_0_10 =0;

data_0 = []

for i in 1:800
    b_spks_0_1, b_ev_0_1 = burst_detection(g_0_1[i],5.0)
    mean_b_spks_0_1 +=b_spks_0_1/800;
    mean_b_ev_0_1 += b_ev_0_1/800;

    b_spks_0_2, b_ev_0_2 = burst_detection(g_0_2[i],5.0)
    mean_b_spks_0_2 +=b_spks_0_2/800;
    mean_b_ev_0_2 += b_ev_0_2/800;

    b_spks_0_3, b_ev_0_3 = burst_detection(g_0_3[i],5.0)
    mean_b_spks_0_3 +=b_spks_0_3/800;
    mean_b_ev_0_3 += b_ev_0_3/800;

    b_spks_0_4, b_ev_0_4 = burst_detection(g_0_4[i],5.0)
    mean_b_spks_0_4 +=b_spks_0_4/800;
    mean_b_ev_0_4 += b_ev_0_4/800;

    b_spks_0_5, b_ev_0_5 = burst_detection(g_0_5[i],5.0)
    mean_b_spks_0_5 +=b_spks_0_5/800;
    mean_b_ev_0_5 += b_ev_0_5/800;

    b_spks_0_6, b_ev_0_6 = burst_detection(g_0_6[i],5.0)
    mean_b_spks_0_6 +=b_spks_0_6/800;
    mean_b_ev_0_6 += b_ev_0_6/800;

    b_spks_0_7, b_ev_0_7 = burst_detection(g_0_7[i],5.0)
    mean_b_spks_0_7 +=b_spks_0_7/800;
    mean_b_ev_0_7 += b_ev_0_7/800;

    b_spks_0_8, b_ev_0_8 = burst_detection(g_0_8[i],5.0)
    mean_b_spks_0_8 +=b_spks_0_8/800;
    mean_b_ev_0_8 += b_ev_0_8/800;

    b_spks_0_9, b_ev_0_9 = burst_detection(g_0_9[i],5.0)
    mean_b_spks_0_9 +=b_spks_0_9/800;
    mean_b_ev_0_9 += b_ev_0_9/800;

    b_spks_0_10, b_ev_0_10 = burst_detection(g_0_10[i],5.0)
    mean_b_spks_0_10 +=b_spks_0_10/800;
    mean_b_ev_0_10 += b_ev_0_10/800;


    push!(data_0,b_ev_0_1)
    push!(data_0,b_ev_0_2)
    push!(data_0,b_ev_0_3)
    push!(data_0,b_ev_0_4)
    push!(data_0,b_ev_0_5)
    push!(data_0,b_ev_0_6)
    push!(data_0,b_ev_0_7)
    push!(data_0,b_ev_0_8)
    push!(data_0,b_ev_0_9)
    push!(data_0,b_ev_0_10)
end

mean_b_spks_0_vector = [mean_b_spks_0_1,mean_b_spks_0_2,mean_b_spks_0_3,mean_b_spks_0_4,mean_b_spks_0_5,mean_b_spks_0_6,mean_b_spks_0_7,mean_b_spks_0_8,mean_b_spks_0_9,mean_b_spks_0_10]
mean_b_ev_0_vector = [mean_b_ev_0_1,mean_b_ev_0_2,mean_b_ev_0_3,mean_b_ev_0_4,mean_b_ev_0_5,mean_b_ev_0_6,mean_b_ev_0_7,mean_b_ev_0_8,mean_b_ev_0_9,mean_b_ev_0_10]
## g=25
steps1 = readdlm("SOM NETWORK/g=25/I=250/steps/Spike_steps_Inject_250_SOM_1.txt",Int64)
index1 = readdlm("SOM NETWORK/g=25/I=250/index/Spike_index_Inject_250_SOM_1.txt",Int64)

steps2 = readdlm("SOM NETWORK/g=25/I=250/steps/Spike_steps_Inject_250_SOM_2.txt",Int64)
index2 = readdlm("SOM NETWORK/g=25/I=250/index/Spike_index_Inject_250_SOM_2.txt",Int64)

steps3 = readdlm("SOM NETWORK/g=25/I=250/steps/Spike_steps_Inject_250_SOM_3.txt",Int64)
index3 = readdlm("SOM NETWORK/g=25/I=250/index/Spike_index_Inject_250_SOM_3.txt",Int64)

steps4 = readdlm("SOM NETWORK/g=25/I=250/steps/Spike_steps_Inject_250_SOM_4.txt",Int64)
index4 = readdlm("SOM NETWORK/g=25/I=250/index/Spike_index_Inject_250_SOM_4.txt",Int64)

steps5 = readdlm("SOM NETWORK/g=25/I=250/steps/Spike_steps_Inject_250_SOM_5.txt",Int64)
index5 = readdlm("SOM NETWORK/g=25/I=250/index/Spike_index_Inject_250_SOM_5.txt",Int64)

steps6 = readdlm("SOM NETWORK/g=25/I=250/steps/Spike_steps_Inject_250_SOM_6.txt",Int64)
index6 = readdlm("SOM NETWORK/g=25/I=250/index/Spike_index_Inject_250_SOM_6.txt",Int64)

steps7 = readdlm("SOM NETWORK/g=25/I=250/steps/Spike_steps_Inject_250_SOM_7.txt",Int64)
index7 = readdlm("SOM NETWORK/g=25/I=250/index/Spike_index_Inject_250_SOM_7.txt",Int64)

steps8 = readdlm("SOM NETWORK/g=25/I=250/steps/Spike_steps_Inject_250_SOM_8.txt",Int64)
index8 = readdlm("SOM NETWORK/g=25/I=250/index/Spike_index_Inject_250_SOM_8.txt",Int64)

steps9 = readdlm("SOM NETWORK/g=25/I=250/steps/Spike_steps_Inject_250_SOM_9.txt",Int64)
index9 = readdlm("SOM NETWORK/g=25/I=250/index/Spike_index_Inject_250_SOM_9.txt",Int64)

steps10 = readdlm("SOM NETWORK/g=25/I=250/steps/Spike_steps_Inject_250_SOM_10.txt",Int64)
index10 = readdlm("SOM NETWORK/g=25/I=250/index/Spike_index_Inject_250_SOM_10.txt",Int64)

g_25_1,g_25_2,g_25_3,g_25_4,g_25_5,g_25_6,g_25_7,g_25_8,g_25_9,g_25_10=spike_events_all(steps1,index1,steps2,index2,steps3,index3,steps4,index4,steps5,index5,steps6,index6,steps7,index7,steps8,index8,steps9,index9,steps10,index10)

mean_b_spks_25_1=0;
mean_b_ev_25_1 =0;

mean_b_spks_25_2=0;
mean_b_ev_25_2 =0;

mean_b_spks_25_3=0;
mean_b_ev_25_3 =0;

mean_b_spks_25_4=0;
mean_b_ev_25_4 =0;

mean_b_spks_25_5=0;
mean_b_ev_25_5 =0;

mean_b_spks_25_6=0;
mean_b_ev_25_6 =0;

mean_b_spks_25_7=0;
mean_b_ev_25_7 =0;

mean_b_spks_25_8=0;
mean_b_ev_25_8 =0;

mean_b_spks_25_9=0;
mean_b_ev_25_9 =0;

mean_b_spks_25_10=0;
mean_b_ev_25_10 =0;

upper_limit_25 = zeros(10)
lower_limit_25= zeros(10)

data_25=[]
for i in 1:800
    b_spks_25_1, b_ev_25_1 = burst_detection(g_25_1[i],5.0)
    mean_b_spks_25_1 +=b_spks_25_1/800;
    mean_b_ev_25_1 += b_ev_25_1/800;

    b_spks_25_2, b_ev_25_2 = burst_detection(g_25_2[i],5.0)
    mean_b_spks_25_2 +=b_spks_25_2/800;
    mean_b_ev_25_2 += b_ev_25_2/800;

    b_spks_25_3, b_ev_25_3 = burst_detection(g_25_3[i],5.0)
    mean_b_spks_25_3 +=b_spks_25_3/800;
    mean_b_ev_25_3 += b_ev_25_3/800;

    b_spks_25_4, b_ev_25_4 = burst_detection(g_25_4[i],5.0)
    mean_b_spks_25_4 +=b_spks_25_4/800;
    mean_b_ev_25_4 += b_ev_25_4/800;

    b_spks_25_5, b_ev_25_5 = burst_detection(g_25_5[i],5.0)
    mean_b_spks_25_5 +=b_spks_25_5/800;
    mean_b_ev_25_5 += b_ev_25_5/800;

    b_spks_25_6, b_ev_25_6 = burst_detection(g_25_6[i],5.0)
    mean_b_spks_25_6 +=b_spks_25_6/800;
    mean_b_ev_25_6 += b_ev_25_6/800;

    b_spks_25_7, b_ev_25_7 = burst_detection(g_25_7[i],5.0)
    mean_b_spks_25_7 +=b_spks_25_7/800;
    mean_b_ev_25_7 += b_ev_25_7/800;

    b_spks_25_8, b_ev_25_8 = burst_detection(g_25_8[i],5.0)
    mean_b_spks_25_8 +=b_spks_25_8/800;
    mean_b_ev_25_8 += b_ev_25_8/800;

    b_spks_25_9, b_ev_25_9 = burst_detection(g_25_9[i],5.0)
    mean_b_spks_25_9 +=b_spks_25_9/800;
    mean_b_ev_25_9 += b_ev_25_9/800;

    b_spks_25_10, b_ev_25_10 = burst_detection(g_25_10[i],5.0)
    mean_b_spks_25_10 +=b_spks_25_10/800;
    mean_b_ev_25_10 += b_ev_25_10/800;

    push!(data_25,b_ev_25_1)
    push!(data_25,b_ev_25_2)
    push!(data_25,b_ev_25_3)
    push!(data_25,b_ev_25_4)
    push!(data_25,b_ev_25_5)
    push!(data_25,b_ev_25_6)
    push!(data_25,b_ev_25_7)
    push!(data_25,b_ev_25_8)
    push!(data_25,b_ev_25_9)
    push!(data_25,b_ev_25_10)
end
mean_b_spks_25_vector = [mean_b_spks_25_1,mean_b_spks_25_2,mean_b_spks_25_3,mean_b_spks_25_4,mean_b_spks_25_5,mean_b_spks_25_6,mean_b_spks_25_7,mean_b_spks_25_8,mean_b_spks_25_9,mean_b_spks_25_10]
mean_b_ev_25_vector = [mean_b_ev_25_1,mean_b_ev_25_2,mean_b_ev_25_3,mean_b_ev_25_4,mean_b_ev_25_5,mean_b_ev_25_6,mean_b_ev_25_7,mean_b_ev_25_8,mean_b_ev_25_9,mean_b_ev_25_10]
## g=50
"""
steps1 = readdlm("SOM NETWORK/g=50/I=250/steps/Spike_steps_Inject_250_SOM_1.txt",Int64)
index1 = readdlm("SOM NETWORK/g=50/I=250/index/Spike_index_Inject_250_SOM_1.txt",Int64)

steps2 = readdlm("SOM NETWORK/g=50/I=250/steps/Spike_steps_Inject_250_SOM_2.txt",Int64)
index2 = readdlm("SOM NETWORK/g=50/I=250/index/Spike_index_Inject_250_SOM_2.txt",Int64)

steps3 = readdlm("SOM NETWORK/g=50/I=250/steps/Spike_steps_Inject_250_SOM_3.txt",Int64)
index3 = readdlm("SOM NETWORK/g=50/I=250/index/Spike_index_Inject_250_SOM_3.txt",Int64)

steps4 = readdlm("SOM NETWORK/g=50/I=250/steps/Spike_steps_Inject_250_SOM_4.txt",Int64)
index4 = readdlm("SOM NETWORK/g=50/I=250/index/Spike_index_Inject_250_SOM_4.txt",Int64)

steps5 = readdlm("SOM NETWORK/g=50/I=250/steps/Spike_steps_Inject_250_SOM_5.txt",Int64)
index5 = readdlm("SOM NETWORK/g=50/I=250/index/Spike_index_Inject_250_SOM_5.txt",Int64)

steps6 = readdlm("SOM NETWORK/g=50/I=250/steps/Spike_steps_Inject_250_SOM_6.txt",Int64)
index6 = readdlm("SOM NETWORK/g=50/I=250/index/Spike_index_Inject_250_SOM_6.txt",Int64)

steps7 = readdlm("SOM NETWORK/g=50/I=250/steps/Spike_steps_Inject_250_SOM_7.txt",Int64)
index7 = readdlm("SOM NETWORK/g=50/I=250/index/Spike_index_Inject_250_SOM_7.txt",Int64)

steps8 = readdlm("SOM NETWORK/g=50/I=250/steps/Spike_steps_Inject_250_SOM_8.txt",Int64)
index8 = readdlm("SOM NETWORK/g=50/I=250/index/Spike_index_Inject_250_SOM_8.txt",Int64)

steps9 = readdlm("SOM NETWORK/g=50/I=250/steps/Spike_steps_Inject_250_SOM_9.txt",Int64)
index9 = readdlm("SOM NETWORK/g=50/I=250/index/Spike_index_Inject_250_SOM_9.txt",Int64)

steps10 = readdlm("SOM NETWORK/g=50/I=250/steps/Spike_steps_Inject_250_SOM_10.txt",Int64)
index10 = readdlm("SOM NETWORK/g=50/I=250/index/Spike_index_Inject_250_SOM_10.txt",Int64)

g_50_1,g_50_2,g_50_3,g_50_4,g_50_5,g_50_6,g_50_7,g_50_8,g_50_9,g_50_10=spike_events_all(steps1,index1,steps2,index2,steps3,index3,steps4,index4,steps5,index5,steps6,index6,steps7,index7,steps8,index8,steps9,index9,steps10,index10)
"""
## g=75
steps1 = readdlm("SOM NETWORK/g=75/I=250/steps/Spike_steps_Inject_250_SOM_1.txt",Int64)
index1 = readdlm("SOM NETWORK/g=75/I=250/index/Spike_index_Inject_250_SOM_1.txt",Int64)

steps2 = readdlm("SOM NETWORK/g=75/I=250/steps/Spike_steps_Inject_250_SOM_2.txt",Int64)
index2 = readdlm("SOM NETWORK/g=75/I=250/index/Spike_index_Inject_250_SOM_2.txt",Int64)

steps3 = readdlm("SOM NETWORK/g=75/I=250/steps/Spike_steps_Inject_250_SOM_3.txt",Int64)
index3 = readdlm("SOM NETWORK/g=75/I=250/index/Spike_index_Inject_250_SOM_3.txt",Int64)

steps4 = readdlm("SOM NETWORK/g=75/I=250/steps/Spike_steps_Inject_250_SOM_4.txt",Int64)
index4 = readdlm("SOM NETWORK/g=75/I=250/index/Spike_index_Inject_250_SOM_4.txt",Int64)

steps5 = readdlm("SOM NETWORK/g=75/I=250/steps/Spike_steps_Inject_250_SOM_5.txt",Int64)
index5 = readdlm("SOM NETWORK/g=75/I=250/index/Spike_index_Inject_250_SOM_5.txt",Int64)

steps6 = readdlm("SOM NETWORK/g=75/I=250/steps/Spike_steps_Inject_250_SOM_6.txt",Int64)
index6 = readdlm("SOM NETWORK/g=75/I=250/index/Spike_index_Inject_250_SOM_6.txt",Int64)

steps7 = readdlm("SOM NETWORK/g=75/I=250/steps/Spike_steps_Inject_250_SOM_7.txt",Int64)
index7 = readdlm("SOM NETWORK/g=75/I=250/index/Spike_index_Inject_250_SOM_7.txt",Int64)

steps8 = readdlm("SOM NETWORK/g=75/I=250/steps/Spike_steps_Inject_250_SOM_8.txt",Int64)
index8 = readdlm("SOM NETWORK/g=75/I=250/index/Spike_index_Inject_250_SOM_8.txt",Int64)

steps9 = readdlm("SOM NETWORK/g=75/I=250/steps/Spike_steps_Inject_250_SOM_9.txt",Int64)
index9 = readdlm("SOM NETWORK/g=75/I=250/index/Spike_index_Inject_250_SOM_9.txt",Int64)

steps10 = readdlm("SOM NETWORK/g=75/I=250/steps/Spike_steps_Inject_250_SOM_10.txt",Int64)
index10 = readdlm("SOM NETWORK/g=75/I=250/index/Spike_index_Inject_250_SOM_10.txt",Int64)

g_75_1,g_75_2,g_75_3,g_75_4,g_75_5,g_75_6,g_75_7,g_75_8,g_75_9,g_75_10=spike_events_all(steps1,index1,steps2,index2,steps3,index3,steps4,index4,steps5,index5,steps6,index6,steps7,index7,steps8,index8,steps9,index9,steps10,index10)

mean_b_spks_75_1=0;
mean_b_ev_75_1 =0;

mean_b_spks_75_2=0;
mean_b_ev_75_2 =0;

mean_b_spks_75_3=0;
mean_b_ev_75_3 =0;

mean_b_spks_75_4=0;
mean_b_ev_75_4 =0;

mean_b_spks_75_5=0;
mean_b_ev_75_5 =0;

mean_b_spks_75_6=0;
mean_b_ev_75_6 =0;

mean_b_spks_75_7=0;
mean_b_ev_75_7 =0;

mean_b_spks_75_8=0;
mean_b_ev_75_8 =0;

mean_b_spks_75_9=0;
mean_b_ev_75_9 =0;

mean_b_spks_75_10=0;
mean_b_ev_75_10 =0;

data_75=[]

for i in 1:800
    b_spks_75_1, b_ev_75_1 = burst_detection(g_75_1[i],5.0)
    mean_b_spks_75_1 +=b_spks_75_1/800;
    mean_b_ev_75_1 += b_ev_75_1/800;

    b_spks_75_2, b_ev_75_2 = burst_detection(g_75_2[i],5.0)
    mean_b_spks_75_2 +=b_spks_75_2/800;
    mean_b_ev_75_2 += b_ev_75_2/800;

    b_spks_75_3, b_ev_75_3 = burst_detection(g_75_3[i],5.0)
    mean_b_spks_75_3 +=b_spks_75_3/800;
    mean_b_ev_75_3 += b_ev_75_3/800;

    b_spks_75_4, b_ev_75_4 = burst_detection(g_75_4[i],5.0)
    mean_b_spks_75_4 +=b_spks_75_4/800;
    mean_b_ev_75_4 += b_ev_75_4/800;

    b_spks_75_5, b_ev_75_5 = burst_detection(g_75_5[i],5.0)
    mean_b_spks_75_5 +=b_spks_75_5/800;
    mean_b_ev_75_5 += b_ev_75_5/800;

    b_spks_75_6, b_ev_75_6 = burst_detection(g_75_6[i],5.0)
    mean_b_spks_75_6 +=b_spks_75_6/800;
    mean_b_ev_75_6 += b_ev_75_6/800;

    b_spks_75_7, b_ev_75_7 = burst_detection(g_75_7[i],5.0)
    mean_b_spks_75_7 +=b_spks_75_7/800;
    mean_b_ev_75_7 += b_ev_75_7/800;

    b_spks_75_8, b_ev_75_8 = burst_detection(g_75_8[i],5.0)
    mean_b_spks_75_8 +=b_spks_75_8/800;
    mean_b_ev_75_8 += b_ev_75_8/800;

    b_spks_75_9, b_ev_75_9 = burst_detection(g_75_9[i],5.0)
    mean_b_spks_75_9 +=b_spks_75_9/800;
    mean_b_ev_75_9 += b_ev_75_9/800;

    b_spks_75_10, b_ev_75_10 = burst_detection(g_75_10[i],5.0)
    mean_b_spks_75_10 +=b_spks_75_10/800;
    mean_b_ev_75_10 += b_ev_75_10/800;

    push!(data_75,b_ev_75_1)
    push!(data_75,b_ev_75_2)
    push!(data_75,b_ev_75_3)
    push!(data_75,b_ev_75_4)
    push!(data_75,b_ev_75_5)
    push!(data_75,b_ev_75_6)
    push!(data_75,b_ev_75_7)
    push!(data_75,b_ev_75_8)
    push!(data_75,b_ev_75_9)
    push!(data_75,b_ev_75_10)

end

mean_b_spks_75_vector = [mean_b_spks_75_1,mean_b_spks_75_2,mean_b_spks_75_3,mean_b_spks_75_4,mean_b_spks_75_5,mean_b_spks_75_6,mean_b_spks_75_7,mean_b_spks_75_8,mean_b_spks_75_9,mean_b_spks_75_10]
mean_b_ev_75_vector = [mean_b_ev_75_1,mean_b_ev_75_2,mean_b_ev_75_3,mean_b_ev_75_4,mean_b_ev_75_5,mean_b_ev_75_6,mean_b_ev_75_7,mean_b_ev_75_8,mean_b_ev_75_9,mean_b_ev_75_10]

## g=150

steps1 = readdlm("SOM NETWORK/g=150/I=250/steps/Spike_steps_Inject_250_SOM_1.txt",Int64)
index1 = readdlm("SOM NETWORK/g=150/I=250/index/Spike_index_Inject_250_SOM_1.txt",Int64)

steps2 = readdlm("SOM NETWORK/g=150/I=250/steps/Spike_steps_Inject_250_SOM_2.txt",Int64)
index2 = readdlm("SOM NETWORK/g=150/I=250/index/Spike_index_Inject_250_SOM_2.txt",Int64)

steps3 = readdlm("SOM NETWORK/g=150/I=250/steps/Spike_steps_Inject_250_SOM_3.txt",Int64)
index3 = readdlm("SOM NETWORK/g=150/I=250/index/Spike_index_Inject_250_SOM_3.txt",Int64)

steps4 = readdlm("SOM NETWORK/g=150/I=250/steps/Spike_steps_Inject_250_SOM_4.txt",Int64)
index4 = readdlm("SOM NETWORK/g=150/I=250/index/Spike_index_Inject_250_SOM_4.txt",Int64)

steps5 = readdlm("SOM NETWORK/g=150/I=250/steps/Spike_steps_Inject_250_SOM_5.txt",Int64)
index5 = readdlm("SOM NETWORK/g=150/I=250/index/Spike_index_Inject_250_SOM_5.txt",Int64)

steps6 = readdlm("SOM NETWORK/g=150/I=250/steps/Spike_steps_Inject_250_SOM_6.txt",Int64)
index6 = readdlm("SOM NETWORK/g=150/I=250/index/Spike_index_Inject_250_SOM_6.txt",Int64)

steps7 = readdlm("SOM NETWORK/g=150/I=250/steps/Spike_steps_Inject_250_SOM_7.txt",Int64)
index7 = readdlm("SOM NETWORK/g=150/I=250/index/Spike_index_Inject_250_SOM_7.txt",Int64)

steps8 = readdlm("SOM NETWORK/g=150/I=250/steps/Spike_steps_Inject_250_SOM_8.txt",Int64)
index8 = readdlm("SOM NETWORK/g=150/I=250/index/Spike_index_Inject_250_SOM_8.txt",Int64)

steps9 = readdlm("SOM NETWORK/g=150/I=250/steps/Spike_steps_Inject_250_SOM_9.txt",Int64)
index9 = readdlm("SOM NETWORK/g=150/I=250/index/Spike_index_Inject_250_SOM_9.txt",Int64)

steps10 = readdlm("SOM NETWORK/g=150/I=250/steps/Spike_steps_Inject_250_SOM_10.txt",Int64)
index10 = readdlm("SOM NETWORK/g=150/I=250/index/Spike_index_Inject_250_SOM_10.txt",Int64)

g_150_1,g_150_2,g_150_3,g_150_4,g_150_5,g_150_6,g_150_7,g_150_8,g_150_9,g_150_10=spike_events_all(steps1,index1,steps2,index2,steps3,index3,steps4,index4,steps5,index5,steps6,index6,steps7,index7,steps8,index8,steps9,index9,steps10,index10)

mean_b_spks_150_1=0;
mean_b_ev_150_1 =0;

mean_b_spks_150_2=0;
mean_b_ev_150_2 =0;

mean_b_spks_150_3=0;
mean_b_ev_150_3 =0;

mean_b_spks_150_4=0;
mean_b_ev_150_4 =0;

mean_b_spks_150_5=0;
mean_b_ev_150_5 =0;

mean_b_spks_150_6=0;
mean_b_ev_150_6 =0;

mean_b_spks_150_7=0;
mean_b_ev_150_7 =0;

mean_b_spks_150_8=0;
mean_b_ev_150_8 =0;

mean_b_spks_150_9=0;
mean_b_ev_150_9 =0;

mean_b_spks_150_10=0;
mean_b_ev_150_10 =0;

data_150=[]
for i in 1:800
    b_spks_150_1, b_ev_150_1 = burst_detection(g_150_1[i],5.0)
    mean_b_spks_150_1 +=b_spks_150_1/800;
    mean_b_ev_150_1 += b_ev_150_1/800;

    b_spks_150_2, b_ev_150_2 = burst_detection(g_150_2[i],5.0)
    mean_b_spks_150_2 +=b_spks_150_2/800;
    mean_b_ev_150_2 += b_ev_150_2/800;


    b_spks_150_3, b_ev_150_3 = burst_detection(g_150_3[i],5.0)
    mean_b_spks_150_3 +=b_spks_150_3/800;
    mean_b_ev_150_3 += b_ev_150_3/800;

    b_spks_150_4, b_ev_150_4 = burst_detection(g_150_4[i],5.0)
    mean_b_spks_150_4 +=b_spks_150_4/800;
    mean_b_ev_150_4 += b_ev_150_4/800;

    b_spks_150_5, b_ev_150_5 = burst_detection(g_150_5[i],5.0)
    mean_b_spks_150_5 +=b_spks_150_5/800;
    mean_b_ev_150_5 += b_ev_150_5/800;

    b_spks_150_6, b_ev_150_6 = burst_detection(g_150_6[i],5.0)
    mean_b_spks_150_6 +=b_spks_150_6/800;
    mean_b_ev_150_6 += b_ev_150_6/800;

    b_spks_150_7, b_ev_150_7 = burst_detection(g_150_7[i],5.0)
    mean_b_spks_150_7 +=b_spks_150_7/800;
    mean_b_ev_150_7 += b_ev_150_7/800;

    b_spks_150_8, b_ev_150_8 = burst_detection(g_150_8[i],5.0)
    mean_b_spks_150_8 +=b_spks_150_8/800;
    mean_b_ev_150_8 += b_ev_150_8/800;

    b_spks_150_9, b_ev_150_9 = burst_detection(g_150_9[i],5.0)
    mean_b_spks_150_9 +=b_spks_150_9/800;
    mean_b_ev_150_9 += b_ev_150_9/800;

    b_spks_150_10, b_ev_150_10 = burst_detection(g_150_10[i],5.0)
    mean_b_spks_150_10 +=b_spks_150_10/800;
    mean_b_ev_150_10 += b_ev_150_10/800;

    push!(data_150,b_ev_150_1)
    push!(data_150,b_ev_150_2)
    push!(data_150,b_ev_150_3)
    push!(data_150,b_ev_150_4)
    push!(data_150,b_ev_150_5)
    push!(data_150,b_ev_150_6)
    push!(data_150,b_ev_150_7)
    push!(data_150,b_ev_150_8)
    push!(data_150,b_ev_150_9)
    push!(data_150,b_ev_150_10)
end

mean_b_spks_150_vector = [mean_b_spks_150_1,mean_b_spks_150_2,mean_b_spks_150_3,mean_b_spks_150_4,mean_b_spks_150_5,mean_b_spks_150_6,mean_b_spks_150_7,mean_b_spks_150_8,mean_b_spks_150_9,mean_b_spks_150_10]
mean_b_ev_150_vector = [mean_b_ev_150_1,mean_b_ev_150_2,mean_b_ev_150_3,mean_b_ev_150_4,mean_b_ev_150_5,mean_b_ev_150_6,mean_b_ev_150_7,mean_b_ev_150_8,mean_b_ev_150_9,mean_b_ev_150_10]

## Plot
mean_bursts_0= sum(mean_b_ev_0_vector)/10;
upper_limit_0 = sqrt(sum((data_0.-mean_bursts_0).^2)/length(data_0))
lower_limit_0 = sqrt(sum((data_0.-mean_bursts_0).^2)/length(data_0))

mean_bursts_25= sum(mean_b_ev_25_vector)/10;
upper_limit_25 = sqrt(sum((data_25.-mean_bursts_25).^2)/length(data_25))
lower_limit_25 = sqrt(sum((data_25.-mean_bursts_25).^2)/length(data_25))

mean_bursts_75= sum(mean_b_ev_75_vector)/10;
upper_limit_75 = sqrt(sum((data_75.-mean_bursts_75).^2)/length(data_75))
lower_limit_75 = sqrt(sum((data_75.-mean_bursts_75).^2)/length(data_75))

mean_bursts_150= sum(mean_b_ev_150_vector)/10;
upper_limit_150 = sqrt(sum((data_150.-mean_bursts_150).^2)/length(data_150))
lower_limit_150 = sqrt(sum((data_150.-mean_bursts_150).^2)/length(data_150))

x = [0,25,75,150];
y = [mean_bursts_0,mean_bursts_25,mean_bursts_75,mean_bursts_150];
yerr=[[lower_limit_0,lower_limit_25,lower_limit_75,lower_limit_150],[upper_limit_0,upper_limit_25,upper_limit_75,upper_limit_150]]
my_dpi=96
fig,ax = PyPlot.subplots(figsize=(1200/my_dpi, 800/my_dpi), dpi=my_dpi)
ax.scatter(x,y,c="k")
ax.errorbar(x,y, yerr=yerr,c="k",capsize=8)
ax.set_ylabel("Burst Events",fontsize=25)
ax.set_xlabel("AIS T-type Conductance (nS)",fontsize=25)
ax.spines["top"].set_visible(false) # Hide the top edge of the axis
ax.spines["right"].set_visible(false)
fig.tight_layout()
=#
