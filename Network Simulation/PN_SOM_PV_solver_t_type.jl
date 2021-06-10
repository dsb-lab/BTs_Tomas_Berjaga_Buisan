function PN_SOM_PV_solve_T_type(par_time::Tuple,par_number::NTuple{3,Int64},par_initial::NTuple{4,Float64},g_Ca_t::Float64,Matrix_connectivity::Array{Float64});
    tvec,dt,Dmax,delay_time = par_time;
    number_pyramidal,number_PV,number_SOM = par_number;
    I_inject_S_initial,I_inject_D_initial,I_inject_PV_initial,I_inject_SOM_initial = par_initial;

    time_length = length(tvec);
    sum_1 = number_pyramidal+number_PV;

    steps_delay=delay_time_syn/dt;
    delay_back_propagation= round(3.0/dt);

    mu_PV = I_inject_PV_initial; mu_SOM = I_inject_SOM_initial;
    mu_S=I_inject_S_initial; mu_D=I_inject_D_initial;

    #Pyramidal Definition
    D=zeros(number_pyramidal);
    state_1=zeros(number_pyramidal);
    state_2=zeros(number_pyramidal);
    state_3=zeros(number_pyramidal);
    state_4=zeros(number_pyramidal);
    state_5=zeros(number_pyramidal);
    state_7=zeros(number_pyramidal);
    state_8=zeros(number_pyramidal);
    state_9=zeros(number_pyramidal);
    state_10=zeros(number_pyramidal);
    state_11=zeros(number_pyramidal);

    m = zeros(number_pyramidal);
    h = zeros(number_pyramidal);
    m_t = zeros(number_pyramidal);
    h_t = zeros(number_pyramidal);

    S_AHP  = zeros(number_pyramidal);
    I_Inject_S = zeros(number_pyramidal);
    I_Inject_D = zeros(number_pyramidal);

    V_S = zeros(time_length,number_pyramidal);
    V_D = zeros(time_length,number_pyramidal);

    #PV Definition
    state_1_PV=zeros(number_PV);
    state_2_PV=zeros(number_PV);
    state_3_PV=zeros(number_PV);
    state_4_PV=zeros(number_PV);

    I_Inject_PV = zeros(number_PV);
    v_PV = zeros(number_PV);
    u_PV = zeros(number_PV);

    #SOM Definition
    state_1_SOM=zeros(number_SOM);
    state_2_SOM=zeros(number_SOM);
    state_3_SOM=zeros(number_SOM);
    state_4_SOM=zeros(number_SOM);

    I_Inject_SOM = zeros(number_SOM);
    v_SOM = zeros(number_SOM);
    u_SOM = zeros(number_SOM);

    #Synapses Definition
    S_syn  = zeros(number_total);
    spkDel= Vector{Int}();
    spkneuronindx= Vector{Int}();

    #Solving Definition
    spk_V_D_step= Vector{Int}();
    spk_V_D_index= Vector{Int}();

    spk_step_pyramidal= Vector{Int}();
    spk_index_pyramidal= Vector{Int}();

    spk_step_all_net= Vector{Int}();
    spk_index_all_net= Vector{Int}();

## Inizialization
for i = 1:number_total
    if i <= number_pyramidal
        V_S[1,i] = V_rest_S;
        V_D[1,i]= V_rest_D;
        I_Inject_S[i] = I_inject_S_initial;
        I_Inject_D[i] = I_inject_D_initial;
        m[i]= 1.0 / (1.0+exp(Sm*(V_rest_D-m_med)));
        h[i]= 1.0 / (1.0+exp(Sh*(V_rest_D-h_med)));
        m_t[i]= 1.0 / (1.0+exp(-(V_rest_S-m_med_t)/k1))
        h_t[i]= 1.0 / (1.0+exp((V_rest_S-h_med_t)/k2))
        state_1[i]=V_rest_S;
   elseif i > number_pyramidal && i <= sum_1
       v_PV[i-number_pyramidal]= -55.0; #mV
       I_Inject_PV[i-number_pyramidal]= I_inject_PV_initial;
   else
       v_SOM[i-sum_1]= v_rest_SOM;
       I_Inject_SOM[i-sum_1]= I_inject_SOM_initial;
   end
end
## Solver
    for (i,t) in enumerate(tvec[1:end-1])
       overprint("step = $i out of $time_length", i=i)
       spkDel.-=1;
       idxs= spkDel.==0;
       idxsyn = spkneuronindx[idxs];
       for tres in idxsyn
          if tres <=number_pyramidal
              S_syn[tres]+=g_syne;
           else
              S_syn[tres]+=g_syni;
           end
       end

       deleteat!(spkDel,idxs);
       deleteat!(spkneuronindx,idxs);

       spk_V_D_step.-=1;
       idxs= spk_V_D_step.==0;
       idxs_V_D = spk_V_D_index[idxs];
       V_D[i,idxs_V_D].+=10.0;
       deleteat!(spk_V_D_step,idxs);
       deleteat!(spk_V_D_index,idxs);

        for neuron in 1:number_pyramidal
            state_2[neuron] = V_D[i,neuron];
            state_3[neuron] = m[neuron];
            state_4[neuron] = h[neuron];
            state_5[neuron] = S_AHP[neuron];
            state_7[neuron] = I_Inject_S[neuron];
            state_8[neuron] = I_Inject_D[neuron];
            state_9[neuron] = S_syn[neuron];
            state_10[neuron] = m_t[neuron];
            state_11[neuron] = h_t[neuron];

            I_Ca = g_Ca *state_3[neuron] * state_4[neuron] * (E_Ca-state_2[neuron]);

            I_AHP = state_5[neuron]*(E_K-state_1[neuron]);
            S_AHP[neuron] = state_5[neuron] - dt *((i_tau_k)*state_5[neuron]);

            S_syn_add_pyr_pyr_dend= Matrix_connectivity[neuron,1:number_pyramidal]'*S_syn[1:number_pyramidal];
            S_syn_add_SOM_pyr_dend = Matrix_connectivity[neuron,sum_1:end]'*S_syn[sum_1:end];

            I_add_pyr_dend = S_syn_add_pyr_pyr_dend*(state_2[neuron]-E_syn_act) + S_syn_add_SOM_pyr_dend*(state_2[neuron]-E_syn_inactivation);

            S_syn[neuron]=state_9[neuron]- dt *(i_tau_syne*state_9[neuron]);

            if D[neuron] == 0.0
                state_1[neuron] = V_S[i,neuron];

                I_t = g_Ca_t*(state_10[neuron]^2.0)*state_11[neuron]*(E_Ca_t-state_1[neuron]);

                S_syn_add_PV_pyr_soma= Matrix_connectivity[neuron,number_pyramidal:sum_1]'*S_syn[number_pyramidal:sum_1];
                I_syn_add_soma = S_syn_add_PV_pyr_soma*(state_1[neuron]-E_syn_inactivation);

                dV_S1 = (i_C_S)*((i_R_S)*(V_rest_S-state_1[neuron])+(i_R_T)*(state_2[neuron]-state_1[neuron])+state_7[neuron]+I_AHP-I_syn_add_soma+I_t) ::Float64;
                V_S[i+1,neuron] = state_1[neuron] + dt*dV_S1;

                if V_S[i+1,neuron] >= threshold
                    state_1[neuron] = V_S[i,neuron] ;
                    S_AHP[neuron] += g_AHP;
                    V_S[i+1,neuron] = V_S_peak ;
                    D[neuron] = Dmax;

                    push!(spkDel,steps_delay)
                    push!(spkneuronindx,neuron)

                    push!(spk_V_D_step,delay_back_propagation)
                    push!(spk_V_D_index,neuron)

                    push!(spk_step_pyramidal,i)
                    push!(spk_index_pyramidal,neuron)

                    push!(spk_step_all_net,i)
                    push!(spk_index_all_net,neuron)

                else
                    state_1[neuron] = V_S[i+1,neuron];
                end
            else
                D[neuron] -= 1.0;
                if D[neuron] == 0
                    V_S[i+1,neuron] = -52.0 ;
                else
                    V_S[i+1,neuron] = V_S_peak;
                end
            end

            dV_D1 = (i_C_D)*(i_R_D*(V_rest_D-state_2[neuron])+i_R_T*(state_1[neuron]-state_2[neuron])+state_8[neuron]+I_Ca-I_add_pyr_dend);
            V_D[i+1,neuron] =  V_D[i,neuron] + dt*dV_D1;

            m_infinite = 1.0 / (1.0+exp(Sm*(state_2[neuron]-m_med)));
            h_infinite = 1.0 / (1.0+exp(Sh*(state_2[neuron]-h_med)));
            m[neuron] = state_3[neuron] + dt*((m_infinite - state_3[neuron]) / tau_m);
            h[neuron] = state_4[neuron] + dt*((h_infinite - state_4[neuron]) / tau_h);

            I_Inject_S[neuron] = state_7[neuron] + ((mu_S-state_7[neuron])/ tau)*dt + sigma*rand(Normal(0, 1))*sqrt(i_tau_2*dt);
            I_Inject_D[neuron] = state_8[neuron] + ((mu_D-state_8[neuron])/ tau)*dt + sigma*rand(Normal(0, 1))*sqrt(i_tau_2*dt);


            tau_m_t = 0.612 + 1 / (exp(-(state_1[neuron]+132.0)/16.7) + exp((state_1[neuron]+16.8)/18.2));

            if state_1[neuron]<-70
                 tau_h_t = exp((state_1[neuron]+467.0)/66.6)
            elseif state_1[neuron]>= -70
                 tau_h_t = exp((state_1[neuron]+22.0)/-10.5) + 28.0
            end

            m_infinite_t = 1.0 / (1.0+exp(-(state_1[neuron]-m_med_t)/k1));
            h_infinite_t = 1.0 / (1.0+exp((state_1[neuron]-h_med_t)/k2));
            m_t[neuron] = state_10[neuron] + dt*((m_infinite_t - state_10[neuron]) / tau_m_t);
            h_t[neuron] = state_11[neuron] + dt*((h_infinite_t - state_11[neuron]) / tau_h_t);
        end

        for PV in 1:number_PV
            state_1_PV[PV] = v_PV[PV];
            state_2_PV[PV] = u_PV[PV]; 
            state_3_PV[PV] = I_Inject_PV[PV];
            state_4_PV[PV] = S_syn[PV+number_pyramidal];

            S_syn_add_PV_PV  = Matrix_connectivity[PV+number_pyramidal,number_pyramidal:sum_1]'*S_syn[number_pyramidal:sum_1];
            S_syn_add_SOM_PV = Matrix_connectivity[PV+number_pyramidal,sum_1:end]'*S_syn[sum_1:end];
            S_syn_add_pyr_PV = Matrix_connectivity[PV+number_pyramidal,1:number_pyramidal]'*S_syn[1:number_pyramidal];

            I_syn_add = S_syn_add_pyr_PV*(state_1_PV[PV]-E_syn_act) + (S_syn_add_SOM_PV+S_syn_add_PV_PV)*(state_1_PV[PV]-E_syn_inactivation);

            S_syn[PV+number_pyramidal]=state_4_PV[PV]- dt *((i_tau_syni)*state_4_PV[PV]);

            if state_1_PV[PV]>=v_b
                U = 0.025*((state_1_PV[PV]-v_b)^3.0) ::Float64;
            else
                U = 0.0;
            end

            dv_D1 = (0.05)*((state_1_PV[PV]+55.0)*(state_1_PV[PV]+40.0)-state_2_PV[PV]+state_3_PV[PV]-I_syn_add);
            v_PV[PV] = state_1_PV[PV] + dt*dv_D1;

            if v_PV[PV] >= threshold_PV
                v_PV[PV] = -45.0 ;
                push!(spkDel,steps_delay)
                push!(spkneuronindx,PV+number_pyramidal)

                push!(spk_step_all_net,i)
                push!(spk_index_all_net,PV+number_pyramidal)
            end

            du_D1 = 0.2*(U-state_2_PV[PV]);
            u_PV[PV]= state_2_PV[PV]+dt*du_D1;
            I_Inject_PV[PV] = state_3_PV[PV] + ((mu_PV-state_3_PV[PV])/ tau)*dt + sigma*rand(Normal(0, 1))*sqrt(i_tau_2*dt);

        end

        for SOM in 1:number_SOM
            state_1_SOM[SOM] = v_SOM[SOM];
            state_2_SOM[SOM] = u_SOM[SOM];
            state_3_SOM[SOM] = I_Inject_SOM[SOM];
            state_4_SOM[SOM] = S_syn[sum_1+SOM];
        #add where
            S_syn_add_PV_SOM = Matrix_connectivity[sum_1+SOM,number_pyramidal:sum_1]'*S_syn[number_pyramidal:sum_1];
            S_syn_add_SOM_SOM = Matrix_connectivity[sum_1+SOM, sum_1:end]'*S_syn[sum_1:end];
            S_syn_add_pyr_SOM = Matrix_connectivity[sum_1+SOM,1:number_pyramidal]'*S_syn[1:number_pyramidal];

            I_syn_add = S_syn_add_pyr_SOM*(state_1_SOM[SOM]-E_syn_act) + (S_syn_add_SOM_SOM+S_syn_add_PV_SOM) *(state_1_SOM[SOM]-E_syn_inactivation) ::Float64;

            S_syn[sum_1+SOM]=state_4_SOM[SOM]- dt *((i_tau_syni)*state_4_SOM[SOM]);

            dv_D1 = (0.01)*((state_1_SOM[SOM]+56.0)*(state_1_SOM[SOM]+42.0)-state_2_SOM[SOM]+state_3_SOM[SOM]-I_syn_add);
            v_SOM[SOM] = state_1_SOM[SOM] + dt*dv_D1;

            du_D1 = 0.03*(8.0*(state_1_SOM[SOM]+56.0)-state_2_SOM[SOM]);
            u_SOM[SOM]= state_2[SOM]+dt*du_D1;

            I_Inject_SOM[SOM] = state_3_SOM[SOM] + ((mu_SOM-state_3_SOM[SOM])/ tau)*dt + sigma*rand(Normal(0, 1))*sqrt(i_tau_2*dt);

            if v_SOM[SOM]>= (40.0-0.1*state_2_SOM[SOM])
                v_SOM[SOM]=-53.0+0.04*state_2[SOM];
                u_SOM[SOM]= min(state_2_SOM[SOM]+20.0,670.0);
                push!(spkDel,steps_delay)
                push!(spkneuronindx,SOM+sum_1)

                push!(spk_step_all_net,i)
                push!(spk_index_all_net,SOM+sum_1)
            end
        end
    end
    a = rand(1:800);
    return spk_step_pyramidal, spk_index_pyramidal,spk_step_all_net,spk_index_all_net,V_S[:,a],V_D[:,a]
end
