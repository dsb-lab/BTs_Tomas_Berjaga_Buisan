function spike_events_all(steps1::Array{Int64,2},index1::Array{Int64,2},steps2::Array{Int64,2},index2::Array{Int64,2},steps3::Array{Int64,2},index3::Array{Int64,2},steps4::Array{Int64,2},index4::Array{Int64,2},steps5::Array{Int64,2},index5::Array{Int64,2},steps6::Array{Int64,2},index6::Array{Int64,2},steps7::Array{Int64,2},index7::Array{Int64,2},steps8::Array{Int64,2},index8::Array{Int64,2},steps9::Array{Int64,2},index9::Array{Int64,2},steps10::Array{Int64,2},index10::Array{Int64,2})
    tvec = collect(0.0:0.1:3000.0);

    spike_events1 = [Vector{Float64}(undef,0) for _ = 1:1000];
    spike_events2 = [Vector{Float64}(undef,0) for _ = 1:1000];
    spike_events3 = [Vector{Float64}(undef,0) for _ = 1:1000];
    spike_events4 = [Vector{Float64}(undef,0) for _ = 1:1000];
    spike_events5 = [Vector{Float64}(undef,0) for _ = 1:1000];
    spike_events6 = [Vector{Float64}(undef,0) for _ = 1:1000];
    spike_events7 = [Vector{Float64}(undef,0) for _ = 1:1000];
    spike_events8 = [Vector{Float64}(undef,0) for _ = 1:1000];
    spike_events9 = [Vector{Float64}(undef,0) for _ = 1:1000];
    spike_events10 = [Vector{Float64}(undef,0) for _ = 1:1000];

    for i in 1:800
       idxs1= index1.==i;
       idxneu1 = steps1[idxs1];

       idxs2= index2.==i;
       idxneu2 = steps2[idxs2];

       idxs3= index3.==i;
       idxneu3 = steps3[idxs3];

       idxs4= index4.==i;
       idxneu4 = steps4[idxs4];

       idxs5= index5.==i;
       idxneu5 = steps5[idxs5];

       idxs6= index6.==i;
       idxneu6 = steps6[idxs6];

       idxs7= index7.==i;
       idxneu7 = steps7[idxs7];

       idxs8= index8.==i;
       idxneu8 = steps8[idxs8];

       idxs9= index9.==i;
       idxneu9 = steps9[idxs9];

       idxs10= index10.==i;
       idxneu10 = steps10[idxs10];

       for j in idxneu1
            push!(spike_events1[i],tvec[j])
       end
       for j in idxneu2
            push!(spike_events2[i],tvec[j])
       end
       for j in idxneu3
            push!(spike_events3[i],tvec[j])
       end
       for j in idxneu4
            push!(spike_events4[i],tvec[j])
       end
       for j in idxneu5
            push!(spike_events5[i],tvec[j])
       end
       for j in idxneu6
            push!(spike_events6[i],tvec[j])
       end
       for j in idxneu7
            push!(spike_events7[i],tvec[j])
       end
       for j in idxneu8
            push!(spike_events8[i],tvec[j])
       end
       for j in idxneu9
            push!(spike_events9[i],tvec[j])
       end
       for j in idxneu10
            push!(spike_events10[i],tvec[j])
       end
    end
    return spike_events1,spike_events2,spike_events3,spike_events4,spike_events5,spike_events6,spike_events7,spike_events8,spike_events9,spike_events10
end

function burst_detection(ts::Vector{Float64}, threshold::Float64; only_times=false)
    delete_spikes!(ts)
    ISI    = computeISI(ts)
    if length(ISI)<1
        if only_times
            return Float64[]
        else
            return 0.0, 0.0, Float64[], 0.0, Float64[]
        end
    end
    b      = ISI .< threshold
    b_spks = sum(b)
    idx_b  = findall(b)
    b_ev   = 0
    b_ev_times = Float64[]
    for i=1:length(idx_b)
        j = idx_b[i]
        if j-1 == 0
            b_ev += 1
            b_spks +=1
        elseif ISI[j-1] > threshold
            push!(b_ev_times, ts[j])
            b_ev += 1
            b_spks +=1
        end
    end
    nonbspks = length(ts) - b_spks
    if only_times
        return b_ev_times
    end
    return b_spks, b_ev# b_ev_times, nonbspks, ISI
end

function delete_spikes!(ts)
    last_idx = 0
    for i=1:length(ts)
        if ts[i] > 500.0
            last_idx = i-1
            break
        end
    end
    deleteat!(ts, 1:last_idx)
end


function computeISI(ts::Vector{Float64})
    if length(ts)<1
        return []
    end
    ISI = zeros(length(ts)-1)
    for i in eachindex(ISI)
        ISI[i] = ts[i+1] - ts[i]
    end
    return ISI
end
