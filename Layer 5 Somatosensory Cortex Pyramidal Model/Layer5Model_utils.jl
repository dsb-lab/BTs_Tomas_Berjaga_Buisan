function pyr_solver(tvec,dt,g, I_Inject_D,g_Ca_t)
    V_S= zeros(length(tvec));
    V_D = zeros(length(tvec));
    m = zeros(length(tvec));
    h = zeros(length(tvec));
    S_AHP  = zeros(length(tvec));
    II_AHP = zeros(length(tvec));
    II_Ca = zeros(length(tvec));
    II_Ca_t = zeros(length(tvec));
    I_poisson = zeros(length(tvec));
    D=0
    #S_AHP[1]=5.3857192733681245;

    V_S[1] = V_rest_S;
    V_D[1] = V_rest_D;
    m[1]= 1 / (1+exp(Sm*(V_rest_D-m_med)));
    h[1]= 1 / (1+exp(Sh*(V_rest_D-h_med)));

    m_t = zeros(length(tvec));
    h_t = zeros(length(tvec));
    m_t[1]= 1 / (1+exp(-(V_rest_S-m_med_t)/k1))
    h_t[1]= 1 / (1+exp((V_rest_S-h_med_t)/k2))

    spike_times_soma = [];
    state = [V_S[1], V_D[1], m[1], h[1], S_AHP[1], m_t[1], h_t[1] ];

    for (i,t) in enumerate(tvec[1:end-1])
        state[2] = V_D[i];
        state[3] = m[i];
        state[4] = h[i];
        state[5] = S_AHP[i];
        state[6] = m_t[i];
        state[7] = h_t[i];

        I_Ca = g_Ca*state[3]*state[4]*(E_Ca-state[2]);
        II_Ca[i] = I_Ca;

        I_t = g_Ca_t*(state[6]^2)*state[7]*(E_Ca_t-state[1]);
        II_Ca_t[i] = I_t

        I_AHP = state[5]*(E_K-state[1]);
        II_AHP[i] = I_AHP;
        S_AHP[i+1] = state[5] - dt *((1/tau_k)*state[5]);
        I_poisson[i] = (g[i]*-state[1])

        if D==0
            state[1] = V_S[i];
            dV_S1 = (1/C_S)*((1/R_S)*(V_rest_S-state[1])+(1/R_T)*(state[2]-state[1])+(g[i]*-state[1])+250+I_AHP+I_t);
            V_S[i+1] = state[1] + dt*dV_S1;

            if V_S[i+1] >= threshold
                state[1] = V_S[i];
                S_AHP[i+1] += g_AHP ;
                V_S[i+1] = V_S_peak;
                D = Dmax
                spike_times_soma= append!(spike_times_soma,t)
            else
                state[1] = V_S[i+1] ;
            end
        else
            D -= 1
            if D==0
                V_S[i+1] = -52.0;
            else
                V_S[i+1] = V_S_peak;
            end
        end

        if length(spike_times_soma)!= 0 && i-round(spike_times_soma[end]/dt) == round(3/dt)
                V_D[i]+=10 ;
        end

        dV_D1 = (1/C_D)*((1/R_D)*(V_rest_D-state[2])+(1/R_T)*(state[1]-state[2])+I_Inject_D+I_Ca);
        V_D[i+1] =  V_D[i] + dt*dV_D1;

        m_infinite = 1 / (1+exp(Sm*(state[2]-m_med)));
        h_infinite = 1 / (1+exp(Sh*(state[2]-h_med)));
        m[i+1] = m[i]+ dt*((m_infinite - state[3]) / tau_m);
        h[i+1] = h[i]+ dt*((h_infinite - state[4]) / tau_h);

        tau_m_t = 0.612 + 1 / (exp(-(state[1]+132)/16.7) + exp((state[1]+16.8)/18.2))

        if state[1]<-80
             tau_h_t = exp((state[1]+467)/66.6)
        elseif state[1]>= -80
             tau_h_t = exp((state[1]+22)/-10.5) + 28
        end

        m_infinite_t = 1 / (1+exp(-(state[1]-m_med_t)/k1));
        h_infinite_t = 1 / (1+exp((state[1]-h_med_t)/k2));
        m_t[i+1] = m_t[i]+ dt*((m_infinite_t - state[6]) / tau_m_t);
        h_t[i+1] = h_t[i]+ dt*((h_infinite_t - state[7]) / tau_h_t);
    end
    return V_S,V_D,II_Ca_t,I_poisson,II_Ca,II_AHP
end
