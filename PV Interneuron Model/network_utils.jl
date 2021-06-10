#Network Functions
function connectivity_network(par_number::NTuple{3,Int64},par_connectivity::NTuple{9,Float64})
       pyramidal,PV,SOM = par_number;
       prob_connect_pyr_pyr,prob_connect_pyr_PV,prob_connect_pyr_SOM,prob_connect_PV_PV,prob_connect_PV_pyr, prob_connect_PV_SOM,prob_connect_SOM_SOM,prob_connect_SOM_pyr,prob_connect_SOM_PV=par_connectivity;
       total_neurons = pyramidal+PV+SOM;
       connectivity_matrix=zeros(total_neurons,total_neurons);
       for j = 1:pyramidal
           for i = 1:pyramidal
               if i!=j && rand() <= prob_connect_pyr_pyr
                   connectivity_matrix[i,j]=1.0;
               end
           end
           for i = pyramidal:(PV+pyramidal)
               if i!=j && rand() <= prob_connect_pyr_PV
                   connectivity_matrix[i,j]=1.0;
               end
           end

           for i = (PV+pyramidal):total_neurons
               if i!=j && rand() <=  prob_connect_pyr_SOM
                   connectivity_matrix[i,j]=1.0;
               end
           end
       end

       for j = pyramidal:(pyramidal+PV)
           for i = 1:pyramidal
               if i!=j && rand() <= prob_connect_PV_pyr
                   connectivity_matrix[i,j]=1.0;
               end
           end
           for i = pyramidal:(PV+pyramidal)
               if i!=j && rand() <= prob_connect_PV_PV
                   connectivity_matrix[i,j]=1.0;
               end
           end

           for i = (PV+pyramidal):total_neurons
               if i!=j && rand() <=  prob_connect_PV_SOM
                   connectivity_matrix[i,j]=1.0;
               end
           end
       end

       for j = (pyramidal+PV):total_neurons
           for i = 1:pyramidal
               if i!=j && rand() <= prob_connect_SOM_pyr
                   connectivity_matrix[i,j]=1.0;
               end
           end
           for i = pyramidal:(PV+pyramidal)
               if i!=j && rand() <= prob_connect_SOM_PV
                   connectivity_matrix[i,j]=1.0;
               end
           end

           for i = (PV+pyramidal):total_neurons
               if i!=j && rand() <=  prob_connect_SOM_SOM
                   connectivity_matrix[i,j]=1.0;
               end
           end
       end
       return connectivity_matrix
   end

   function raster_plot(spk_index::Vector{Int},spk_step::Vector{Int},number::Int,tvec::Vector{Float64})
       spike_events = [Vector{Float64}(undef,0) for _ = 1:number]
       for i in 1:number
          idxs= spk_index.==i;
          idxneu = spk_step[idxs];
          for j in idxneu
               push!(spike_events[i],tvec[j])
           end
       end
       return spike_events
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

   function fftfreqs(n, d) #longitud tvec, dt
       if isodd(n)
           left  = range(0, stop=(n-1)/2, step=1)
           right = range(-(n-1)/2, stop=-1,step=1)
           freqs = [left ; right]
       elseif iseven(n)
           left  = range(0, stop=n/2 - 1, step=1)
           right = range(-n/2, stop=-1,step=1)
           freqs = [left ; right]
       end
       return freqs./(d*n)
   end

   function PS(dt::Float64, data::Vector{Float64}; only_pos=true)
       n  = length(data)
       n2 = floor(Int, n/2)
       d  = dt ./1000.0
       PS = abs.(fft(data)).^2
       PS ./= n
       fr = fftfreqs(n, d)
       return PS[2:n2], fr[2:n2]
   end

   function overprint(str::String; i=1)
       if i==1
           println(str)
       elseif i>1
           print("\u1b[1F")
           #Moves cursor to beginning of the line n (default 1) lines up
           print(str)   #prints the new line
           print("\u1b[0K")
           # clears  part of the line.
           #If n is 0 (or missing), clear from cursor to the end of the line.
           #If n is 1, clear from cursor to beginning of the line.
           #If n is 2, clear entire line.
           #Cursor position does not change.
           println() #prints a new line, i really don't like this arcane codes
       end
   end
