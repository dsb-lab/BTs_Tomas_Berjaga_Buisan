function gauss(times::Vector{Float64}, sig::Float64, mn::Float64)
    dist = zeros(length(times))
    for i in eachindex(dist)
        dist[i] += (1/(sig*2pi)) * exp(-0.5 * ((times[i]-mn)/sig)^2)
    end
    return dist
end

function gaussians(simtime::Float64, h::Float64, sig::Float64, freq::Float64; dephase=0.0)
    times      = collect(range(0 , step=h, stop=simtime))
    per        = 1000.0/freq
    cycles     = round(Int, simtime/per,RoundDown)
    centers    = zeros(cycles)
    deph       = per*dephase
    distsum = zeros(length(times))
    for c in eachindex(centers)
        centers[c]    = per*c + deph
        distsum .+= gauss(times, sig, centers[c])
    end
    return distsum
end

function osc_poisson(dt::Float64, simtime::Float64, steps::Int64, params::Vector{Float64}; g=false, dephase=0.0)
    sig    = params[1]
    freq   = params[2]
    gain   = params[3]
    offset = params[4]
    gain /= freq
    dist  = gain .* gaussians(simtime, dt, sig, freq, dephase=dephase) .+ offset
    raster = zeros(Int8, steps)
    samps = rand(steps)
    for i=1:steps
        if samps[i] < dist[i]
            raster[i] = 1
        end
    end
    if g
        return raster, dist
    end
    return raster
end

function PSC(dt::Float64,ts::Vector,tvec::Vector)
   S = zeros(length(tvec))
    for i=2:length(tvec)
        dS = -S[i-1]*tau_syne
        S[i] = S[i-1] + dt*dS
        nspks = sum(ts.== i)
        while nspks > 0
            S[i] += g_syne
            nspks -=1
        end
    end
    return S
end

function compute_ts(raster)
    ts = Vector{Float64}()
    for i in eachindex(raster)
        if raster[i] == 1
            push!(ts, i)
        end
    end
    return ts
end

function fftfreqs(n, d)
    if isodd(n)
        left  = range(0, stop=(n-1)/2, step=1)
        right = range((n-1)/2, stop=-1,step=1)
        freqs = [left ; right]
    elseif iseven(n)
        left  = range(0, stop=n/2 - 1, step=1)
        right = range(-n/2, stop=-1,step=1)
        freqs = [left ; right]
    end
    return freqs./(d*n)
end

function PS(h::Float64, data::Vector; only_pos=true)
    n  = length(data)
    n2 = floor(Int, n/2)
    d  = h ./1000.0
    PS = abs.(fft(data)).^2
    PS ./= n
    fr = fftfreqs(n, d)
    return PS[2:n2], fr[2:n2]
end
