#Global parameters pyramidal
    const R_T= 0.065; #ohms
    const i_R_T = 1.0/R_T;

    const R_S= 0.05;  #ohms
    const i_R_S = 1.0/R_S;

    const R_D= 0.043;  #ohms
    const i_R_D = 1.0/R_S;

    const C_S= 260.0;  #pF
    const i_C_S = 1.0/C_S;

    const C_D= 120.0;   #pF
    const i_C_D = 1.0/C_D;

    const V_rest_S=-70.0 ; #mV
    const V_rest_D= -60.0;   #mV
    const V_S_peak = 10.0;
    const threshold= -47.0;  #mV

    const g_AHP=4.0;  #nS
    const g_Ca=70.0;  #nS
    const E_K=-90.0;  #mV
    const Sm= -0.5;  #mV^-1
    const Sh= 0.5;  #mV^-1
    const m_med=-9.0; #mV
    const h_med=-21.0;  #mV
    const tau_m=15.0;  #ms
    const tau_h=80.0;  #ms

    const tau_k=80.0;  #ms
    const i_tau_k = 1.0/tau_k;

    const E_Ca = 120.0; #mV

#Global parameters PV interneurons
    const v_b = -55.0;  #mV
    const threshold_PV = 25.0;  #mV

    const step_current_PV=300.0; #mV
    const init_step_function_time=50.0;  #ms
    const final_step_function_time=125.0;  #ms

#Global parameters SOM interneurons
    const v_rest_SOM=-56.0;  #mV

#Global parameters for synapses
    const tau_syne= 2.5;
    const i_tau_syne = 1.0/tau_syne;

    const tau_syni= 5.0;
    const i_tau_syni = 1.0/tau_syni;

    const g_syne=2.0;
    const g_syni=3.0;
    const E_syn_act=0.0;
    const E_syn_inactivation=-75.0;
    const delay_time_syn=2.0; #ms

#Global Parameters T-type Channels
#Ornhstein-Uhlenbeck parameters
    const sigma=300.0; const tau =3.0;
    const i_tau_2 = 2.0/tau;
