using Distributions
using Random
using DataFrames
using DataFramesMeta
using CSV
using Plots

###################
#### FUNCTIONS ####
###################

function p_exp(I_t, N_t, β)
  p_exp = 1-exp(β * - (I_t/N_t))
  return(p_exp)
end

#################
#### PROCESS ####
#################

Random.seed!(4)

function SEIRsim(S0, E0, I0, R0, T, β, γ, δ)

  data = fill(-1., T, 5)
  track = fill(-1., T, 4)

  St = S0; Et = E0; It = I0; Rt = R0
  Nt = St + Et + It + Rt

  p_inf = 1 - exp(-γ)
  p_rem = 1 - exp(-δ)

  for t in 1:T

    p_exp_t = p_exp(It, Nt, β)

    new_E = rand(Binomial(St, p_exp_t))
    new_I = rand(Binomial(Et, p_inf))
    new_R = rand(Binomial(It, p_rem))

    St = St - new_E
    Et = Et + new_E - new_I
    It = It + new_I - new_R
    Rt = Rt + new_R
    Nt = St + Et + It + Rt

    Pt = rand(Binomial(20, (It/Nt)))

    data[t, :] = [St, Et, It, Rt, Pt]
    track[t, :] = [new_E, new_I, new_R, p_exp_t]

  end

  return(data, track)
end


#################################
#### GENERATE SIMULATED DATA ####
#################################

S0=999; E0=0; I0=1; R0=0; T=200

β_tr=0.25; γ_tr=0.08; δ_tr=0.22

Random.seed!(40)
data_sim_core, track_sim_core = SEIRsim(S0, E0, I0, R0, T, β_tr, γ_tr, δ_tr)

plot(data_sim_core[:, 1:4], label = ["S" "E" "I" "R"], legend = :left, color = [:green :orange :red :blue])

CSV.write("Data/simulated_data_1000.csv", DataFrame(data_sim_core, [:S, :E, :I, :R, :ProportionI]), header = false)
CSV.write("Data/simulated_track_1000.csv", DataFrame(track_sim_core, [:NewE, :NewI, :NewR, :p_exp_t]), header = false)

#################################
#### TEST FOR MODE AROUND 25 ####
#################################

sum_I = fill(-10, 100000)
seeds = fill(-10, 100000)
for i in 1:100000
  Random.seed!(i)
  data_sim_core, track_sim_core = SEIRsim(S0, E0, I0, R0, T, β_tr, γ_tr, δ_tr)
  sum_I[i] = data_sim_core[T, 4]
  if (data_sim_core[T, 4] == 25)
    seeds[i] = i
  end
end
seeds = seeds[seeds .> 0]
println(seeds[1:10])
histogram(sum_I)
