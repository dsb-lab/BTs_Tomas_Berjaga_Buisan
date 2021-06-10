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

#T-type
    const E_Ca_t = 120.0;  #Mirar
    const k1 = 5.2;
    const k2 = 6.2;
    const m_med_t = -48.4; #CaV 3.2
    const h_med_t = -75.6;
