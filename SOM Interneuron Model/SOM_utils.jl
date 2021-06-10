function step_current(dt::Float64,tvec::Vector{Float64},initial::Int,final::Int,Amplitude::Int)
    I = zeros(length(tvec));
    initial = convert(Int,round(initial/dt)); final =convert(Int,round(final/dt));
    for j in initial:final
        I[j]=Amplitude;
    end
    return I
end

function simulate_SOM(tvec::Vector{Float64},v::Vector{Float64},u::Vector{Float64},I::Vector{Float64})
    spike_times=zeros(length(tvec))
    for (i,t) in enumerate(tvec[1:end-1])
        state_1 = v[i]
        state_2 = u[i]
        dv_D1 = (1/100)*((state_1+56)*(state_1+42)-state_2+I[i]);
        v[i+1] = state_1 + dt*dv_D1

        du_D1 = 0.03*(8*(state_1+56)-state_2);
        u[i+1]= state_2+dt*du_D1;

        if v[i+1]>= (40-0.1*state_2)
            v[i+1]=-53+0.04*state_2;
            u[i+1]= min(state_2+20,670);
            spike_times[i+1]=1
        end
    end
    return v,u,I,spike_times
end

function nullcline(v::Vector{Float64},I,u)
    #First Nullcline v=0
    u_nullcline = zeros(length(v));
    u_nullcline_2 = zeros(length(u));
    u_nullcline_3 = zeros(length(u));
    for (i,t) in enumerate(v)
        if v[i]>=v_b
            U = 0.025*((v[i]-v_b)^3)
        else
            U = 0;
        end
        u_nullcline[i] = U+I
    end

    for (i,t) in enumerate(u)
        u_nullcline_2[i] = 0.5 * (-95-sqrt(4*u[i]+225-4*I))
        u_nullcline_3[i] = 0.5 * (-95+sqrt(4*u[i]+225-4*I))
    end
    return u_nullcline,u_nullcline_2,u_nullcline_3
end


function frequency(dt::Float64, tvec::Vector{Float64},v::Vector{Float64},u::Vector{Float64},I::Vector{Int64})
    frequency = []
    length_time = length(I)
    for i in I
        overprint("step = $i out of $length_time", i=i)
        I = step_current(dt::Float64,tvec::Vector{Float64},300::Int,1700::Int,i)
        v,u,I,spikes=simulate_SOM(tvec::Vector{Float64},v::Vector{Float64},u::Vector{Float64},I::Vector{Float64})
        initial = convert(Int,round(500/dt)); final =convert(Int,round(1500/dt));
        frequency_1 = sum(spikes[initial:final])
        push!(frequency,frequency_1)
    end
    return frequency
end
