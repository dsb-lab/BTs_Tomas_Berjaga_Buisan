function gating_variables(m_med::Float64,h_med::Float64,k1::Float64,k2::Float64,voltage::Vector{Float64})
    m=zeros(length(voltage));
    h=zeros(length(voltage));
    for (i,t) in enumerate(voltage)
        m[i] = (1 / (1+exp(-(t-m_med)/k1)))^2;
        h[i] = 1 / (1+exp((t-h_med)/k2));
    end
    return m,h
end

function compute_tau_m(voltage::Vector{Float64})
    tau_m=zeros(length(voltage));
    for (i,t) in enumerate(voltage)
        tau_m[i] = 0.612 + 1 / (exp(-(t+132)/16.7) + exp((t+16.8)/18.2))
    end
    return tau_m
end

function compute_tau_h(voltage::Vector{Float64})
    tau_h=zeros(length(voltage));
    for (i,t) in enumerate(voltage)
        if t<-70
             tau_h[i] = exp((t+467)/66.6)
        elseif t>= -70
             tau_h[i] = exp((t+22)/-10.5) + 28
        end
    end
    return tau_h
end
