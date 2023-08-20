using Distributions
using Random
using Plots
using DataFrames
using DataFramesMeta
using CSV
using JLD2

################
#### SET UP ####
################

include("SEIR - MCMC.jl")

# data_sim_core = Array(CSV.read("Data/simulated_data_1000.csv", DataFrame, header = false))
# track_sim_core = Array(CSV.read("Data/simulated_track_1000.csv", DataFrame, header = false))

init_conds_core = Array(CSV.read("Data/init_conds_agg_point2.csv", DataFrame, header = false))
data_sim_core = Array(CSV.read("Data/simulated_data_1000_agg_point2.csv", DataFrame, header = false))
track_sim_core = Array(CSV.read("Data/simulated_track_1000_agg_point2.csv", DataFrame, header = false))

init_conds_core = Array(CSV.read("Data/init_conds_agg_1.csv", DataFrame, header = false))
data_sim_core = Array(CSV.read("Data/simulated_data_1000_agg_1.csv", DataFrame, header = false))
track_sim_core = Array(CSV.read("Data/simulated_track_1000_agg_1.csv", DataFrame, header = false))

init_conds_core = Array(CSV.read("Data/init_conds_agg_7.csv", DataFrame, header = false))
data_sim_core = Array(CSV.read("Data/simulated_data_1000_agg_7.csv", DataFrame, header = false))
track_sim_core = Array(CSV.read("Data/simulated_track_1000_agg_7.csv", DataFrame, header = false))

init_conds_core = Array(CSV.read("Data/init_conds_agg_30.csv", DataFrame, header = false))
data_sim_core = Array(CSV.read("Data/simulated_data_1000_agg_30.csv", DataFrame, header = false))
track_sim_core = Array(CSV.read("Data/simulated_track_1000_agg_30.csv", DataFrame, header = false))



plot(data_sim_core[:, 1:4], label = ["S" "E" "I" "R"], legend = :left, color = [:green :orange :red :blue])

# S0=999; E0=0; I0=1; R0=0; T=200
S0=init_conds_core[1]; E0=init_conds_core[2]; I0=init_conds_core[3]; R0=init_conds_core[4]; T=size(data_sim_core,1)
β_tr=0.25; γ_tr=0.08; δ_tr=0.22

d_β_core = Gamma(5, 0.05)
d_γ_core = Gamma(1.6, 0.05)
d_δ_core = Gamma(4.4, 0.05)

# Half mean
# d_β_core = Gamma(2.5, 0.05)
# d_γ_core = Gamma(0.8, 0.05)
# d_δ_core = Gamma(2.2, 0.05)

# d_β_core = Uniform(0, 0.6)
# d_γ_core = Uniform(0, 0.2)
# d_δ_core = Uniform(0, 0.4)

initial_conditions_core = [S0, E0, I0, R0]




###################
#### EXAMPLE 1 ####
###################

data_sim_1 = deepcopy(data_sim_core)
track_sim_1 = deepcopy(track_sim_core)

params_dists_1 = deepcopy([d_β_core, d_γ_core, d_δ_core])
initial_conditions_1 = deepcopy(initial_conditions_core)

# init_values_1 = [β_tr, γ_tr, δ_tr]
init_values_1 = [missing, missing, missing]
data_aug_inf_1 = [false, false, false, false, false, false]

λ_init = 0.3
m_init = 1.4

Random.seed!(3)

@time begin
  res_1, other_res_1, aug_res_1,
   mSE_track_1, mEI_track_1,
   arSE_track_1, arEI_track_1, SIM_track_1, EVENTS_track_1 =
                    Blk_Adaptive_MCMC_post_all(N_its = 2005000, data_init = data_sim_1, track_init = track_sim_1,
                                                prior_dists = params_dists_1, params_true = init_values_1,
                                                λ_init = λ_init, m_init = m_init,
                                                data_aug_infer = data_aug_inf_1, param_infer = true,
                                                init_cond = initial_conditions_1, T = T,
                                                nt_up = [1,1])
end

describe(res_1)
llh(data_sim_1, track_sim_1, [β_tr, γ_tr, δ_tr], initial_conditions_1)

llh(data_sim_1, track_sim_1, [0.0892039, 0.246152, 0.0803677], initial_conditions_1)

plot(Array(res_1[:, 1]))
plot(Array(res_1[:, 2]))
plot(Array(res_1[:, 3]))

# save("Inference/Set 1/params/SIM_track.jld2", "SIM", SIM_track_1)
# save("Inference/Set 1/params/EVENTS_track.jld2", "EVENTS", EVENTS_track_1)
CSV.write("Inference/Set 1/params/res.csv", res_1, header=true)
CSV.write("Inference/Set 1/params/other_res.csv", other_res_1, header=true)





###################
#### EXAMPLE 2 ####
###################

data_sim_2 = deepcopy(data_sim_core)
track_sim_2 = deepcopy(track_sim_core)

params_dists_2 = deepcopy([d_β_core, d_γ_core, d_δ_core])
initial_conditions_2 = deepcopy(initial_conditions_core)

init_values_2 = [missing, missing, missing]
data_aug_inf_2 = [false, false, true, false, false, false]

λ_init = 0.3
m_init = 1.4

Random.seed!(3)

@time begin
  res_2, other_res_2, aug_res_2,
   mSE_track_2, mEI_track_2, mIR_track_2,
   arSE_track_2, arEI_track_2, arIR_track_2, SIM_track_2, EVENTS_track_2 =
                    Blk_Adaptive_MCMC_post_all(N_its = 1005000, data_init = data_sim_2, track_init = track_sim_2,
                                                prior_dists = params_dists_2, params_true = init_values_2,
                                                λ_init = λ_init, m_init = m_init,
                                                data_aug_infer = data_aug_inf_2, param_infer = true,
                                                init_cond = initial_conditions_2, T = T,
                                                nt_up = [1,1])
end


# save("Inference/Set 1/params/SIM_track.jld2", "SIM", SIM_track_2)
# save("Inference/Set 1/params/EVENTS_track.jld2", "EVENTS", EVENTS_track_2)
CSV.write("Inference/Set 1/mSE/res.csv", res_2, header=true)
CSV.write("Inference/Set 1/mSE/other_res.csv", other_res_2, header=true)
CSV.write("Inference/Set 1/mSE/aug_res.csv", aug_res_2, header=true)

CSV.write("Inference/Set 1/mSE/mSE_track.csv", mSE_track_2, header=true)




###################
#### EXAMPLE 3 ####
###################

data_sim_3 = deepcopy(data_sim_core)
track_sim_3 = deepcopy(track_sim_core)

d_β_3 = deepcopy(d_β_core)
d_γ_3 = deepcopy(d_γ_core)
d_δ_3 = deepcopy(d_δ_core)

params_dists_3 = [d_β_3, d_γ_3, d_δ_3]
initial_conditions_3 = deepcopy(initial_conditions_core)

init_values_3 = [missing, missing, missing]
data_aug_inf_3 = [false, false, false, true, false, false]

λ_init = 0.3
m_init = 1.4

Random.seed!(3)

@time begin
  res_3, other_res_3, aug_res_3,
   mSE_track_3, mEI_track_3, mIR_track_3,
   arSE_track_3, arEI_track_3, arIR_track_3, SIM_track_3, EVENTS_track_3 =
                    Blk_Adaptive_MCMC_post_all(N_its = 1005000, data_init = data_sim_3, track_init = track_sim_3,
                                                prior_dists = params_dists_3, params_true = init_values_3,
                                                λ_init = λ_init, m_init = m_init,
                                                data_aug_infer = data_aug_inf_3, param_infer = true,
                                                init_cond = initial_conditions_3, T = T,
                                                nt_up = [1,1])
end

# save("Inference/Set 1/mEI/SIM_track.jld2", "SIM", SIM_track_3)
# save("Inference/Set 1/mEI/EVENTS_track.jld2", "EVENTS", EVENTS_track_3)
CSV.write("Inference/Set 1/mEI/res.csv", res_3, header=true)
CSV.write("Inference/Set 1/mEI/other_res.csv", other_res_3, header=true)
CSV.write("Inference/Set 1/mEI/aug_res.csv", aug_res_3, header=true)

CSV.write("Inference/Set 1/mEI/mEI_track.csv", mEI_track_3, header=true)




###################
#### EXAMPLE 4 ####
###################

data_sim_4 = deepcopy(data_sim_core)
track_sim_4 = deepcopy(track_sim_core)

params_dists_4 = deepcopy([d_β_core, d_γ_core, d_δ_core])
initial_conditions_4 = deepcopy(initial_conditions_core)

init_values_4 = [missing, missing, missing]
data_aug_inf_4 = [false, false, false, false, true, false]

λ_init = 0.3
m_init = 1.4

Random.seed!(3)

@time begin
  res_4, other_res_4, aug_res_4,
   mSE_track_4, mEI_track_4, mIR_track_4,
   arSE_track_4, arEI_track_4, arIR_track_4, SIM_track_4, EVENTS_track_4 =
                    Blk_Adaptive_MCMC_post_all(N_its = 1005000, data_init = data_sim_4, track_init = track_sim_4,
                                                prior_dists = params_dists_4, params_true = init_values_4,
                                                λ_init = λ_init, m_init = m_init,
                                                data_aug_infer = data_aug_inf_4, param_infer = true,
                                                init_cond = initial_conditions_4, T = T,
                                                nt_up = [1,1])
end

# save("Inference/Set 1/params/SIM_track.jld2", "SIM", SIM_track_4)
# save("Inference/Set 1/params/EVENTS_track.jld2", "EVENTS", EVENTS_track_4)
CSV.write("Inference/Set 1/arSE/res.csv", res_4, header=true)
CSV.write("Inference/Set 1/arSE/other_res.csv", other_res_4, header=true)
CSV.write("Inference/Set 1/arSE/aug_res.csv", aug_res_4, header=true)

CSV.write("Inference/Set 1/arSE/arSE_track.csv", arSE_track_4, header=true)





###################
#### EXAMPLE 5 ####
###################

data_sim_5 = deepcopy(data_sim_core)
track_sim_5 = deepcopy(track_sim_core)

params_dists_5 = deepcopy([d_β_core, d_γ_core, d_δ_core])
initial_conditions_5 = deepcopy(initial_conditions_core)

init_values_5 = [missing, missing, missing]
data_aug_inf_5 = [false, false, false, false, false, true]

λ_init = 0.3
m_init = 1.4

Random.seed!(3)

@time begin
  res_5, other_res_5, aug_res_5,
   mSE_track_5, mEI_track_5, mIR_track_5,
   arSE_track_5, arEI_track_5, arIR_track_5, SIM_track_5, EVENTS_track_5 =
                    Blk_Adaptive_MCMC_post_all(N_its = 1005000, data_init = data_sim_5, track_init = track_sim_5,
                                                prior_dists = params_dists_5, params_true = init_values_5,
                                                λ_init = λ_init, m_init = m_init,
                                                data_aug_infer = data_aug_inf_5, param_infer = true,
                                                init_cond = initial_conditions_5, T = T,
                                                nt_up = [1,1])
end


# save("Inference/Set 1/params/SIM_track.jld2", "SIM", SIM_track_5)
# save("Inference/Set 1/params/EVENTS_track.jld2", "EVENTS", EVENTS_track_5)
CSV.write("Inference/Set 1/arEI/res.csv", res_5, header=true)
CSV.write("Inference/Set 1/arEI/other_res.csv", other_res_5, header=true)
CSV.write("Inference/Set 1/arEI/aug_res.csv", aug_res_5, header=true)

CSV.write("Inference/Set 1/arEI/arEI_track.csv", arEI_track_5, header=true)






###################
#### EXAMPLE 6 ####
###################

data_sim_6 = deepcopy(data_sim_core)
track_sim_6 = deepcopy(track_sim_core)

d_β_6 = deepcopy(d_β_core)
d_γ_6 = deepcopy(d_γ_core)
d_δ_6 = deepcopy(d_δ_core)

params_dists_6 = [d_β_6, d_γ_6, d_δ_6]
initial_conditions_6 = deepcopy(initial_conditions_core)

init_values_6 = [missing, missing, missing]
data_aug_inf_6 = [false, false, true, true, false, false]

λ_init = 0.3
m_init = 1.4

Random.seed!(3)

@time begin
  res_6, other_res_6, aug_res_6,
   mSE_track_6, mEI_track_6, mIR_track_6,
   arSE_track_6, arEI_track_6, arIR_track_6, SIM_track_6, EVENTS_track_6 =
                    Blk_Adaptive_MCMC_post_all(N_its = 1005000, data_init = data_sim_6, track_init = track_sim_6,
                                                prior_dists = params_dists_6, params_true = init_values_6,
                                                λ_init = λ_init, m_init = m_init,
                                                data_aug_infer = data_aug_inf_6, param_infer = true,
                                                init_cond = initial_conditions_6, T = T,
                                                nt_up = [1,1])
end

# save("Inference/Set 1/mSEmEI/SIM_track.jld2", "SIM", SIM_track_6)
# save("Inference/Set 1/mSEmEI/EVENTS_track.jld2", "EVENTS", EVENTS_track_6)
CSV.write("Inference/Set 1/mSEmEI/res.csv", res_6, header=true)
CSV.write("Inference/Set 1/mSEmEI/other_res.csv", other_res_6, header=true)
CSV.write("Inference/Set 1/mSEmEI/aug_res.csv", aug_res_6, header=true)

CSV.write("Inference/Set 1/mSEmEI/mSE_track.csv", mSE_track_6, header=true)
CSV.write("Inference/Set 1/mSEmEI/mEI_track.csv", mEI_track_6, header=true)



###################
#### EXAMPLE 7 ####
###################

data_sim_7 = deepcopy(data_sim_core)
track_sim_7 = deepcopy(track_sim_core)

d_β_7 = deepcopy(d_β_core)
d_γ_7 = deepcopy(d_γ_core)
d_δ_7 = deepcopy(d_δ_core)

params_dists_7 = [d_β_7, d_γ_7, d_δ_7]
initial_conditions_7 = deepcopy(initial_conditions_core)

init_values_7 = [missing, missing, missing]
data_aug_inf_7 = [false, false, false, false, true, true]

λ_init = 0.3
m_init = 1.4

Random.seed!(3)

@time begin
  res_7, other_res_7, aug_res_7,
   mSE_track_7, mEI_track_7, mIR_track_7,
   arSE_track_7, arEI_track_7, arIR_track_7, SIM_track_7, EVENTS_track_7 =
                    Blk_Adaptive_MCMC_post_all(N_its = 1005000, data_init = data_sim_7, track_init = track_sim_7,
                                                prior_dists = params_dists_7, params_true = init_values_7,
                                                λ_init = λ_init, m_init = m_init,
                                                data_aug_infer = data_aug_inf_7, param_infer = true,
                                                init_cond = initial_conditions_7, T = T,
                                                nt_up = [1,1])
end

# save"Inference/Set 1/arSEarEI/SIM_track.jld2", "SIM", SIM_track_7)
# save("Inference/Set 1/arSEarEI/EVENTS_track.jld2", "EVENTS", EVENTS_track_7)
CSV.write("Inference/Set 1/arSEarEI/res.csv", res_7, header=true)
CSV.write("Inference/Set 1/arSEarEI/other_res.csv", other_res_7, header=true)
CSV.write("Inference/Set 1/arSEarEI/aug_res.csv", aug_res_7, header=true)

CSV.write("Inference/Set 1/arSEarEI/arSE_track.csv", arSE_track_7, header=true)
CSV.write("Inference/Set 1/arSEarEI/arEI_track.csv", arEI_track_7, header=true)




###################
#### EXAMPLE 8 ####
###################

data_sim_8 = deepcopy(data_sim_core)
track_sim_8 = deepcopy(track_sim_core)

d_β_8 = deepcopy(d_β_core)
d_γ_8 = deepcopy(d_γ_core)
d_δ_8 = deepcopy(d_δ_core)

params_dists_8 = [d_β_8, d_γ_8, d_δ_8]
initial_conditions_8 = deepcopy(initial_conditions_core)

init_values_8 = [missing, missing, missing]
data_aug_inf_8 = [false, false, true, false, true, false]

λ_init = 0.3
m_init = 1.4

Random.seed!(3)

@time begin
  res_8, other_res_8, aug_res_8,
   mSE_track_8, mEI_track_8, mIR_track_8,
   arSE_track_8, arEI_track_8, arIR_track_8, SIM_track_8, EVENTS_track_8 =
                    Blk_Adaptive_MCMC_post_all(N_its = 1005000, data_init = data_sim_8, track_init = track_sim_8,
                                                prior_dists = params_dists_8, params_true = init_values_8,
                                                λ_init = λ_init, m_init = m_init,
                                                data_aug_infer = data_aug_inf_8, param_infer = true,
                                                init_cond = initial_conditions_8, T = T,
                                                nt_up = [1,1])
end


# save("Inference/Set 1/mSEarSE/SIM_track.jld2", "SIM", SIM_track_8)
# save("Inference/Set 1/mSEarSE/EVENTS_track.jld2", "EVENTS", EVENTS_track_8)
CSV.write("Inference/Set 1/mSEarSE/res.csv", res_8, header=true)
CSV.write("Inference/Set 1/mSEarSE/other_res.csv", other_res_8, header=true)
CSV.write("Inference/Set 1/mSEarSE/aug_res.csv", aug_res_8, header=true)

CSV.write("Inference/Set 1/mSEarSE/mSE_track.csv", mSE_track_8, header=true)
CSV.write("Inference/Set 1/mSEarSE/arSE_track.csv", arSE_track_8, header=true)




###################
#### EXAMPLE 9 ####
###################

data_sim_9 = deepcopy(data_sim_core)
track_sim_9 = deepcopy(track_sim_core)

d_β_9 = deepcopy(d_β_core)
d_γ_9 = deepcopy(d_γ_core)
d_δ_9 = deepcopy(d_δ_core)

params_dists_9 = [d_β_9, d_γ_9, d_δ_9]
initial_conditions_9 = deepcopy(initial_conditions_core)

init_values_9 = [missing, missing, missing]
data_aug_inf_9 = [false, false, false, true, false, true]

λ_init = 0.3
m_init = 1.4

Random.seed!(3)

@time begin
  res_9, other_res_9, aug_res_9,
   mSE_track_9, mEI_track_9, mIR_track_9,
   arSE_track_9, arEI_track_9, arIR_track_9, SIM_track_9, EVENTS_track_9 =
                    Blk_Adaptive_MCMC_post_all(N_its = 1005000, data_init = data_sim_9, track_init = track_sim_9,
                                                prior_dists = params_dists_9, params_true = init_values_9,
                                                λ_init = λ_init, m_init = m_init,
                                                data_aug_infer = data_aug_inf_9, param_infer = true,
                                                init_cond = initial_conditions_9, T = T,
                                                nt_up = [1,1])
end


# save("Inference/Set 1/mEIarEI/SIM_track.jld2", "SIM", SIM_track_9)
# save("Inference/Set 1/mEIarEI/EVENTS_track.jld2", "EVENTS", EVENTS_track_9)
CSV.write("Inference/Set 1/mEIarEI/res.csv", res_9, header=true)
CSV.write("Inference/Set 1/mEIarEI/other_res.csv", other_res_9, header=true)
CSV.write("Inference/Set 1/mEIarEI/aug_res.csv", aug_res_9, header=true)

CSV.write("Inference/Set 1/mEIarEI/mSE_track.csv", mSE_track_9, header=true)
CSV.write("Inference/Set 1/mEIarEI/mEI_track.csv", mEI_track_9, header=true)
CSV.write("Inference/Set 1/mEIarEI/arSE_track.csv", arSE_track_9, header=true)
CSV.write("Inference/Set 1/mEIarEI/arEI_track.csv", arEI_track_9, header=true)




####################
#### EXAMPLE 10 ####
####################

data_sim_10 = deepcopy(data_sim_core)
track_sim_10 = deepcopy(track_sim_core)

d_β_10 = deepcopy(d_β_core)
d_γ_10 = deepcopy(d_γ_core)
d_δ_10 = deepcopy(d_δ_core)

params_dists_10 = [d_β_10, d_γ_10, d_δ_10]
initial_conditions_10 = deepcopy(initial_conditions_core)

init_values_10 = [missing, missing, missing]
data_aug_inf_10 = [false, false, true, true, true, true]

λ_init = 0.3
m_init = 1.4

Random.seed!(3)

@time begin
  res_10, other_res_10, aug_res_10,
   mSE_track_10, mEI_track_10, mIR_track_10,
   arSE_track_10, arEI_track_10, arIR_track_10, SIM_track_10, EVENTS_track_10 =
                    Blk_Adaptive_MCMC_post_all(N_its = 1005000, data_init = data_sim_10, track_init = track_sim_10,
                                                prior_dists = params_dists_10, params_true = init_values_10,
                                                λ_init = λ_init, m_init = m_init,
                                                data_aug_infer = data_aug_inf_10, param_infer = true,
                                                init_cond = initial_conditions_10, T = T,
                                                nt_up = [1,1])
end

# save("Inference/Set 1/mSEmEIarSEarEI/SIM_track.jld2", "SIM", SIM_track_10)
# save("Inference/Set 1/mSEmEIarSEarEI/EVENTS_track.jld2", "EVENTS", EVENTS_track_10)
CSV.write("Inference/Set 1/mSEmEIarSEarEI/res.csv", res_10, header=true)
CSV.write("Inference/Set 1/mSEmEIarSEarEI/other_res.csv", other_res_10, header=true)
CSV.write("Inference/Set 1/mSEmEIarSEarEI/aug_res.csv", aug_res_10, header=true)

CSV.write("Inference/Set 1/mSEmEIarSEarEI/mSE_track.csv", mSE_track_10, header=true)
CSV.write("Inference/Set 1/mSEmEIarSEarEI/mEI_track.csv", mEI_track_10, header=true)
CSV.write("Inference/Set 1/mSEmEIarSEarEI/arSE_track.csv", arSE_track_10, header=true)
CSV.write("Inference/Set 1/mSEmEIarSEarEI/arEI_track.csv", arEI_track_10, header=true)


####################
#### EXAMPLE 11 ####
####################

data_sim_11 = deepcopy(data_sim_core)
track_sim_11 = deepcopy(track_sim_core)

params_dists_11 = deepcopy([d_β_core, d_γ_core, d_δ_core])
initial_conditions_11 = deepcopy(initial_conditions_core)

init_values_11 = [missing, missing, missing]
data_aug_inf_11 = [true, true, true, true, true, true]

λ_init = 0.3
m_init = 0.3

Random.seed!(3)

@time begin
  res_11, other_res_11, aug_res_11,
   mSE_track_11, mEI_track_11,
   arSE_track_11, arEI_track_11, SIM_track_11, EVENTS_track_11 =
                    Blk_Adaptive_MCMC_post_all(N_its = 3005000, data_init = data_sim_11, track_init = track_sim_11,
                                                prior_dists = params_dists_11, params_true = init_values_11,
                                                λ_init = λ_init, m_init = m_init,
                                                data_aug_infer = data_aug_inf_11, param_infer = true,
                                                init_cond = initial_conditions_11, T = T,
                                                nt_up = [1,1], timestep_size = 0.2)
end

describe(res_11)

plot(Array(res_11[:, 1]))
plot(Array(res_11[:, 2]))
plot(Array(res_11[:, 3]))

save("Inference/Set 1/agg_point2/SIM_track.jld2", "SIM", SIM_track_11)
save("Inference/Set 1/agg_point2/EVENTS_track.jld2", "EVENTS", EVENTS_track_11)
CSV.write("Inference/Set 1/agg_point2/res.csv", res_11, header=true)
CSV.write("Inference/Set 1/agg_point2/other_res.csv", other_res_11, header=true)
CSV.write("Inference/Set 1/agg_point2/aug_res.csv", aug_res_11, header=true)

CSV.write("Inference/Set 1/agg_point2/mSE_track.csv", mSE_track_11, header=true)
CSV.write("Inference/Set 1/agg_point2/mEI_track.csv", mEI_track_11, header=true)
CSV.write("Inference/Set 1/agg_point2/arSE_track.csv", arSE_track_11, header=true)
CSV.write("Inference/Set 1/agg_point2/arEI_track.csv", arEI_track_11, header=true)




####################
#### EXAMPLE 12 ####
####################

data_sim_12 = deepcopy(data_sim_core)
track_sim_12 = deepcopy(track_sim_core)

params_dists_12 = deepcopy([d_β_core, d_γ_core, d_δ_core])
initial_conditions_12 = deepcopy(initial_conditions_core)

init_values_12 = [missing, missing, missing]
data_aug_inf_12 = [true, true, true, true, true, true]

λ_init = 0.15
m_init = 0.3

Random.seed!(3)

@time begin
  res_12, other_res_12, aug_res_12,
   mSE_track_12, mEI_track_12,
   arSE_track_12, arEI_track_12, SIM_track_12, EVENTS_track_12 =
                    Blk_Adaptive_MCMC_post_all(N_its = 3005000, data_init = data_sim_12, track_init = track_sim_12,
                                                prior_dists = params_dists_12, params_true = init_values_12,
                                                λ_init = λ_init, m_init = m_init,
                                                data_aug_infer = data_aug_inf_12, param_infer = true,
                                                init_cond = initial_conditions_12, T = T,
                                                nt_up = [1,1], timestep_size = 1)
end

describe(res_12)

plot(Array(res_12[:, 1]))
plot(Array(res_12[:, 2]))
plot(Array(res_12[:, 3]))

save("Inference/Set 1/agg_1/SIM_track.jld2", "SIM", SIM_track_12)
save("Inference/Set 1/agg_1/EVENTS_track.jld2", "EVENTS", EVENTS_track_12)
CSV.write("Inference/Set 1/agg_1/res.csv", res_12, header=true)
CSV.write("Inference/Set 1/agg_1/other_res.csv", other_res_12, header=true)
CSV.write("Inference/Set 1/agg_1/aug_res.csv", aug_res_12, header=true)

CSV.write("Inference/Set 1/agg_1/mSE_track.csv", mSE_track_12, header=true)
CSV.write("Inference/Set 1/agg_1/mEI_track.csv", mEI_track_12, header=true)
CSV.write("Inference/Set 1/agg_1/arSE_track.csv", arSE_track_12, header=true)
CSV.write("Inference/Set 1/agg_1/arEI_track.csv", arEI_track_12, header=true)



####################
#### EXAMPLE 13 ####
####################

data_sim_13 = deepcopy(data_sim_core)
track_sim_13 = deepcopy(track_sim_core)

params_dists_13 = deepcopy([d_β_core, d_γ_core, d_δ_core])
initial_conditions_13 = deepcopy(initial_conditions_core)

init_values_13 = [missing, missing, missing]
data_aug_inf_13 = [true, true, true, true, true, true]

λ_init = 0.15
m_init = 0.3

Random.seed!(3)

@time begin
  res_13, other_res_13, aug_res_13,
   mSE_track_13, mEI_track_13,
   arSE_track_13, arEI_track_13, SIM_track_13, EVENTS_track_13 =
                    Blk_Adaptive_MCMC_post_all(N_its = 3005000, data_init = data_sim_13, track_init = track_sim_13,
                                                prior_dists = params_dists_13, params_true = init_values_13,
                                                λ_init = λ_init, m_init = m_init,
                                                data_aug_infer = data_aug_inf_13, param_infer = true,
                                                init_cond = initial_conditions_13, T = T,
                                                nt_up = [1,1], timestep_size = 7)
end

describe(res_13)

plot(Array(res_13[:, 1]))
plot(Array(res_13[:, 2]))
plot(Array(res_13[:, 3]))

save("Inference/Set 1/agg_7/SIM_track.jld2", "SIM", SIM_track_13)
save("Inference/Set 1/agg_7/EVENTS_track.jld2", "EVENTS", EVENTS_track_13)
CSV.write("Inference/Set 1/agg_7/res.csv", res_13, header=true)
CSV.write("Inference/Set 1/agg_7/other_res.csv", other_res_13, header=true)
CSV.write("Inference/Set 1/agg_7/aug_res.csv", aug_res_13, header=true)

CSV.write("Inference/Set 1/agg_7/mSE_track.csv", mSE_track_13, header=true)
CSV.write("Inference/Set 1/agg_7/mEI_track.csv", mEI_track_13, header=true)
CSV.write("Inference/Set 1/agg_7/arSE_track.csv", arSE_track_13, header=true)
CSV.write("Inference/Set 1/agg_7/arEI_track.csv", arEI_track_13, header=true)




####################
#### EXAMPLE 14 ####
####################

data_sim_14 = deepcopy(data_sim_core)
track_sim_14 = deepcopy(track_sim_core)

params_dists_14 = deepcopy([d_β_core, d_γ_core, d_δ_core])
initial_conditions_14 = deepcopy(initial_conditions_core)

init_values_14 = [missing, missing, missing]
data_aug_inf_14 = [true, true, true, true, true, true]

λ_init = 0.15
m_init = 0.3

Random.seed!(3)

@time begin
  res_14, other_res_14, aug_res_14,
   mSE_track_14, mEI_track_14,
   arSE_track_14, arEI_track_14, SIM_track_14, EVENTS_track_14 =
                    Blk_Adaptive_MCMC_post_all(N_its = 3005000, data_init = data_sim_14, track_init = track_sim_14,
                                                prior_dists = params_dists_14, params_true = init_values_14,
                                                λ_init = λ_init, m_init = m_init,
                                                data_aug_infer = data_aug_inf_14, param_infer = true,
                                                init_cond = initial_conditions_14, T = T,
                                                nt_up = [1,1], timestep_size = 30)
end

describe(res_14)

plot(Array(res_14[:, 1]))
plot(Array(res_14[:, 2]))
plot(Array(res_14[:, 3]))

save("Inference/Set 1/agg_30/SIM_track.jld2", "SIM", SIM_track_14)
save("Inference/Set 1/agg_30/EVENTS_track.jld2", "EVENTS", EVENTS_track_14)
CSV.write("Inference/Set 1/agg_30/res.csv", res_14, header=true)
CSV.write("Inference/Set 1/agg_30/other_res.csv", other_res_14, header=true)
CSV.write("Inference/Set 1/agg_30/aug_res.csv", aug_res_14, header=true)

CSV.write("Inference/Set 1/agg_30/mSE_track.csv", mSE_track_14, header=true)
CSV.write("Inference/Set 1/agg_30/mEI_track.csv", mEI_track_14, header=true)
CSV.write("Inference/Set 1/agg_30/arSE_track.csv", arSE_track_14, header=true)
CSV.write("Inference/Set 1/agg_30/arEI_track.csv", arEI_track_14, header=true)
