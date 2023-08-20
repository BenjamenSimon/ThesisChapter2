using Distributions
using InvertedIndices
using Random
using Plots
using DataFrames
using DataFramesMeta
using LinearAlgebra
using CSV
using OnlineStats
using StatsBase

##############################
### LIKELIHOOD + POSTERIOR ###
##############################

function p_exp(I_t, N_t, β, timestep_size)
  p_exp = 1-exp(- β * (I_t/N_t) * timestep_size)
  return(p_exp)
end

function init_cond_llh(init_conditions, params, p_inf, p_rem, track_in, timestep_size)

  p_exp_init = p_exp(init_conditions[3], sum(init_conditions), params[1], timestep_size)

  prob_new_E_1 = logpdf(Binomial(init_conditions[1], p_exp_init), track_in[1, 1])
  prob_new_I_1 = logpdf(Binomial(init_conditions[2], p_inf), track_in[1, 2])
  prob_new_R_1 = logpdf(Binomial(init_conditions[3], p_rem), track_in[1, 3])

  return(prob_new_E_1 + prob_new_I_1 + prob_new_R_1)
end

function bulk_llh(t, data_in, track_in, params, p_inf, p_rem, timestep_size)

  p_exp_t = p_exp(data_in[(t-1), 3], sum(data_in[(t-1), :]), params[1], timestep_size)
  # println("I = ", data_in[(t-1), 3], " and N = ", sum(data_in[(t-1), :]), " and β = ", params[1])
  # println("p_exp_t = ", p_exp_t, " at time ", t)

  prob_new_E = logpdf(Binomial(convert(Int64, data_in[(t-1), 1]), p_exp_t), track_in[t, 1])
  prob_new_I = logpdf(Binomial(convert(Int64, data_in[(t-1), 2]), p_inf), track_in[t, 2])
  # println("E = ", data_in[(t-1), 2], " and newEI = ", track_in[t, 2])
  # println("p_inf = ", p_inf, " at time ", t)
  prob_new_R = logpdf(Binomial(convert(Int64, data_in[(t-1), 3]), p_rem), track_in[t, 3])

  return(prob_new_E + prob_new_I + prob_new_R)
end

function llh(data_in, track_in, params, init_conditions, timestep_size)

  p_inf = 1 - exp(-params[2] * timestep_size)
  p_rem = 1 - exp(-params[3] * timestep_size)

  llh = init_cond_llh(init_conditions, params, p_inf, p_rem, track_in, timestep_size)


  for t in 2:size(data_in, 1)
    llh = llh + bulk_llh(t, data_in, track_in, params, p_inf, p_rem, timestep_size)
  end

  return(llh)
end

function posterior(data_in, track_in, params, prior_dists, init_conditions, timestep_size)

  loglh = llh(data_in, track_in, params, init_conditions, timestep_size)

  priors = logpdf(prior_dists[1], params[1]) +
           logpdf(prior_dists[2], params[2]) +
           logpdf(prior_dists[3], params[3])

  post = loglh + priors

  return(post)
end


# llh(data_sim_core, track_sim_core, [0.25, 0.08, 0.22], [999,0,1,0])


#################
### PROPOSALS ###
#################

function accepted_prop(n_tune, other_res)

  ## Calculate acceptance rate of last batch ##
  batch_start = ((n_tune - 1) * 25) + 1
  batch_end = (n_tune*25) - 1

  # Accepted Proportion
  return( sum(other_res[batch_start:batch_end, 3]) / 24 )
end

function calculate_delta_n(n_tune, acc_prop)

  if acc_prop < 0.33
    delta_n = -1 * min(0.05, 1/sqrt(n_tune))
  else
    delta_n = min(0.05, 1/sqrt(n_tune))
  end

  return(delta_n)
end

function tune_λ(λ_cur, n_tune, other_res)

    acc_prop = accepted_prop(n_tune, other_res)

    ## Update tuning parameter ##
    delta_n = calculate_delta_n(n_tune, acc_prop)

    # log(λ_new) = log(λ_old) + Δn
    # λ_new = exp( log(λ_old) + Δn )
    # λ_new = exp(log(λ_old)) * exp(Δn)
  return(λ_cur * exp(delta_n))
end

#~~~~~#

function mixture_less(d, λ, log_params_cur)

  log_params_draw = rand(MvNormal(log_params_cur, ( (1/d) * (λ)^2 * I ) ))
end

function mixture_more(it, results, covarM, m, log_params_cur)

  log_params_draw = rand(MvNormal(log_params_cur, ((m)^2 * cov(covarM)) ))

  return(log_params_draw)
end

function propose_params(N_its, results, other_res,
                        it, log_params_cur,
                        n_tune, m, λ, covarM)

  ##############################################
  ### Initialise the parameters and epidemic ###
  ##############################################

      # log_params_draw = deepcopy(log_params_cur)

  ######################
  ### Initial tuning ###
  ######################

      if it < min(5000, (N_its/10))

        if it == (n_tune*25)
          λ = tune_λ(λ, n_tune, other_res)
        end

        log_params_draw = mixture_less(3, λ, log_params_cur)

        mixture = 0
      end

  ##########################
  ### Draw the paramters ###
  ##########################

      if it >= min(5000, (N_its/10))

        mixture = rand()

        ### Non-adaptive draw ###

        if mixture < 0.05
          log_params_draw = mixture_less(3, λ, log_params_cur)
        end

        ### Adaptive Draw ###

        if mixture >= 0.05
          log_params_draw = mixture_more(it, results, covarM, m, log_params_cur)
        end

      end

      #log_params_draw, log_q_ratio, mixture, λ
  return(log_params_draw, 0, mixture, λ)
end

#~~~~~#

function propose_init(data_in, track_in, event, params_, T, init_cond_cur, nt_up, timestep_size)

  data_cur = deepcopy(data_in)
  track_cur = deepcopy(track_in)

  data_prime = deepcopy(data_in)
  track_prime = deepcopy(track_in)

  init_conds_prime = deepcopy(init_cond_cur)

  log_q_ratio = 0.

  #### Idenitfy possible changes ####

  possible_changes = Int[]

  # If there is an E in init conds to move
  if init_cond_cur[(event+1)] > 0
    push!(possible_changes, 0)
  end

  # If there is an SE event at t=1 to move
  if track_cur[1, (event)] > 0
    push!(possible_changes, 1)
  end

  # If no possible moves
  if size(possible_changes, 1) == 0
    Move_track = [0, 0, 3, -1, -1]
    return(data_cur, track_cur, init_cond_cur, -Inf, Move_track)
          # data_in, track_in, log_q_ratio, Move_track
  end

  #### Generate the updates ####

  update_direction = rand(possible_changes)


  if update_direction == 0
    init_conds_prime[event] += 1
    init_conds_prime[event+1] += -1

    track_prime[1, event] += 1
  end

  if update_direction == 1
    init_conds_prime[event] += -1
    init_conds_prime[event+1] += 1

    track_prime[1, event] += -1
  end

  ##### Calculate the log q ratio (proposal density ratio) ####

  q_prime_given_cur = 1 / size(possible_changes, 1)


  possible_changes_prime = Int[]

  # If there is an E in init conds to move
  if init_conds_prime[(event+1)] > 0
    push!(possible_changes_prime, 0)
  end

  # If there is an SE event at t=1 to move
  if track_prime[1, (event)] > 0
    push!(possible_changes_prime, 1)
  end

  q_cur_given_prime = 1 / size(possible_changes_prime, 1)

  log_q_ratio += (q_cur_given_prime - q_prime_given_cur)




  Move_track = [0, 0, 0, -1, -1]

    # data_in, track_in, log_q_ratio, Move_track
  return(data_prime, track_prime, init_conds_prime, log_q_ratio, Move_track)
end

#~~~~~#

function move_find_t_valid(track, event, mint, maxt)

  choosing_a_t = Int64[]
  weights_t = Int64[]

  for t in mint:maxt
    if track[t, event] > 0
      push!(choosing_a_t, t)
      push!(weights_t, track[t, event])
    end
  end

  return(choosing_a_t, weights_t)
end

function move_update_states(data, t, event, Δ, num_event_move_t)

  sgnΔ = sign(Δ)
  A = convert(Int64, (1 - sign(Δ))/2)
  B = 1-A

  for k in (t+(A*Δ)):(t-1+(B*Δ))

    data[k , event] += sgnΔ * num_event_move_t

    data[k , (event+1)] -= sgnΔ * num_event_move_t
  end

  return(data)
end

function move_update_track(track, t, event, Δ, num_event_move_t)

  track[t , event] -= num_event_move_t
  track[(t+Δ) , event] += num_event_move_t

  return(track)
end

function move_calc_q_prime_given_cur(t, choosing_a_t, weights_t)

  prob = (weights_t[choosing_a_t .== t] / sum(weights_t))[1]

  # prob = 1/size(choosing_a_t, 1)

  return(log(prob))
end

function move_calc_q_cur_given_prime(t, Δ, track_prime, event, mint, maxt)

  choosing_a_t_prime, weights_t_prime = move_find_t_valid(track_prime, event, mint, maxt)

  prob = (weights_t_prime[choosing_a_t_prime .== (t+Δ)] / sum(weights_t_prime))[1]

  # prob = 1/size(choosing_a_t_prime, 1)

  return(log(prob))
end

function propose_move(data_in, track_in, event, params_, T, init_cond_cur, nt_up, timestep_size)

  data_cur = deepcopy(data_in)
  track_cur = deepcopy(track_in)

  data_prime = deepcopy(data_in)
  track_prime = deepcopy(track_in)

  log_q_ratio = 0.

  #### Generate the updates ####

  Δ_array = rand([-1, 1], nt_up)

  for Δi in Δ_array

    #### Choose a timestep ####

      if Δi == 1 # move forward in time
        choosing_a_t, weights_t = move_find_t_valid(track_cur, event, 1, (T-1))
      end

      if Δi == -1 # move backward in time
        choosing_a_t, weights_t = move_find_t_valid(track_cur, event, 2, T)
      end

    #### Generate update ####

      t = sample(choosing_a_t, weights(weights_t))[1]

      # Number of [] to [] events at time t
      num_event_t = track_cur[t, event]

      # Generate a number to be moved
      num_event_move_t = 1

    #### Generate new states ####

      # Update the events
      track_prime = move_update_track(track_cur, t, event, Δi, num_event_move_t)

      # Update the states
      data_prime = move_update_states(data_cur, t, event, Δi, num_event_move_t)


      #### Quick check for validity ####

      A = convert(Int64, (1 - sign(Δi))/2)
      B = 1-A

      # Early return: Invalid Update
      for j in 1:4
        for k in (t+(A*Δi)):(t+(B*Δi))

          if data_prime[k,j] < 0
            Move_track = [t, 0, 3, Δi, num_event_move_t]
            # data_in, track_in, log_q_ratio, Move_track
            return(data_in, track_in, init_cond_cur, -Inf, Move_track)
          end

        end
      end

      ##### Calculate the log q ratio (proposal density ratio) ####

      q_prime_given_cur = move_calc_q_prime_given_cur(t, choosing_a_t, weights_t)

      if Δi == 1 # move forward in time
        q_cur_given_prime = move_calc_q_cur_given_prime(t, Δi, track_prime, event, 2, T)
      end

      if Δi == -1 # move backward in time
        q_cur_given_prime = move_calc_q_cur_given_prime(t, Δi, track_prime, event, 1, (T-1))
      end

      log_q_ratio += (q_cur_given_prime - q_prime_given_cur)

      #### Update intermediate arrays ####

      data_cur = deepcopy(data_prime)
      track_cur = deepcopy(track_prime)

    end # end of for Δi in Δ_array

  # ELSE

  Move_track = [0, 0, 0, -1, -1]

    # data_in, track_in, log_q_ratio, Move_track
  return(data_prime, track_prime, init_cond_cur, log_q_ratio, Move_track)
end

#~~~~~#

function add_find_t_valid(data, event, T, init_conditions, params_, timestep_size)

  choosing_a_t = Int64[]

  for t in 2:T
    if (data[(t-1), event] > 0. && addrem_calc_probs(data, params_, event, t, init_conditions, timestep_size) > 0.)
      push!(choosing_a_t, t)
    end
  end

  return(choosing_a_t)
end

function rem_find_t_valid(track, event, T)

  choosing_a_t = Int64[]

  for t in 2:T
    if track[t, event] > 0.
      push!(choosing_a_t, t)
    end
  end

  return(choosing_a_t)
end

function addrem_calc_probs(data, params_, event, t, init_conditions, timestep_size)

  if event == 1
    if t == 1
      I_t = init_conditions[3]
      N_t = sum(init_conditions)
    else
      I_t = data[(t-1), 3]
      N_t = sum(data[(t-1), 1:4])
    end

    # p_exp_t
    prob = p_exp(I_t, N_t, params_[1], timestep_size)
  end

  if event == 2
    # p_inf
    prob = 1 - exp(-params_[2])
  end

  if event == 3
    # p_rem
    prob = 1 - exp(-params_[3])
  end

  return(prob)
end

function add_calc_weights(choosing_a_t, data, params_, event, init_conditions, timestep_size)

  wghts = Float64[]

  for t in choosing_a_t
      push!(wghts, addrem_calc_probs(data, params_, event, t, init_conditions, timestep_size))
  end

  return(wghts)
end

function addrem_update_states(data, t, T, event, Δ)

  for k in t:T
    data[k, event] -= Δ
    data[k, (event+1)] += Δ
  end

  return(data)
end

function addrem_update_track(track, t, event, Δ)

  track[t , event] += Δ

  return(track)
end

function add_calc_q_prime_given_cur(t, choosing_a_t, wghts)

  pos = (choosing_a_t .== t)

  prob = wghts[pos][1]/sum(wghts)

  return(log(prob))
end

function add_calc_q_cur_given_prime(t, data_prime, event, T, init_conditions, params_, timestep_size)

    choosing_a_t_prime = add_find_t_valid(data_prime, event, T, init_conditions, params_, timestep_size)

    wghts_prime = add_calc_weights(choosing_a_t_prime, data_prime, params_, event, init_conditions, timestep_size)

    pos = (choosing_a_t_prime .== t)

    prob = wghts_prime[pos][1]/sum(wghts_prime)

    return(log(prob))
end

function rem_calc_q_prime_given_cur(choosing_a_t)

  prob = 1/size(choosing_a_t, 1)

  return(log(prob))
end

function rem_calc_q_cur_given_prime(track_prime, event, T)

  choosing_a_t_prime = rem_find_t_valid(track_prime, event, T)

  prob = 1/size(choosing_a_t_prime, 1)

  return(log(prob))
end

function propose_AddRem(data_in, track_in, event, params_, T, init_cond_cur, nt_up, timestep_size)

  data_cur = deepcopy(data_in)
  track_cur = deepcopy(track_in)

  data_prime = deepcopy(data_in)
  track_prime = deepcopy(track_in)

  log_q_ratio = 0.

  #### Choose to add or remove events

  Δ_array = rand([-1, 1], nt_up)

  #### Generate updates ####

  for Δi in Δ_array

    #### Choose a timestep ####

      if Δi == 1 # add

        choosing_a_t = add_find_t_valid(data_cur, event, T, init_cond_cur, params_, timestep_size)

        if size(choosing_a_t, 1) > 0
          wghts = add_calc_weights(choosing_a_t, data_cur, params_, event, init_cond_cur, timestep_size)
          t = sample(choosing_a_t, Weights(wghts))[1]
        else
          #skip
          AddRem_track = [-1, 0, 3, -1, -1, -1, -1, -1]
          # data_in, track_in, log_q_ratio, Move_track
          return(data_in, track_in, init_cond_cur, -Inf, AddRem_track)
        end
      end

      if Δi == -1 # remove

        choosing_a_t = rem_find_t_valid(track_cur, event, T)

        if size(choosing_a_t, 1) > 0
          t = sample(choosing_a_t, 1, replace = false)[1]
        else
          AddRem_track = [-1, 0, 3, -1, -1, -1, -1, -1]
          # data_in, track_in, log_q_ratio, Move_track
          return(data_in, track_in, init_cond_cur, -Inf, AddRem_track)
        end
      end

    #### Generate update ####

      # Generate a new number of events

      num_event_t = track_in[t, event] # Number of [] to [] events at time t

      new_event_t = num_event_t + Δi

    #### Generate new states ####

      # Update the events
      track_prime = addrem_update_track(track_prime, t, event, Δi)

      # Update the states
      data_prime = addrem_update_states(data_prime, t, T, event, Δi)

      # Early return: Invalid Update
      for j in 1:4
        for k in t:T
          if data_prime[k,j] < 0
            AddRem_track = [t, 0, 3, num_event_t, new_event_t, Δi, -1, -1]
            # data_in, track_in, log_q_ratio, Move_track
            return(data_in, track_in, init_cond_cur, -Inf, AddRem_track)
          end
        end
      end
      # I wonder could we make this a skip and put the early return version
      # after all the updates have been proposed?


    ##### Calculate the log q ratio (proposal density ratio) ####

    # println("current choosing_a_t = ", choosing_a_t)

      if Δi == 1 # add
        q_prime_given_cur = add_calc_q_prime_given_cur(t, choosing_a_t, wghts)
        q_cur_given_prime = rem_calc_q_cur_given_prime(track_prime, event, T)
      end

      if Δi == -1 # remove
        q_prime_given_cur = rem_calc_q_prime_given_cur(choosing_a_t)
        q_cur_given_prime = add_calc_q_cur_given_prime(t, data_prime, event, T, init_cond_cur, params_, timestep_size)
      end

      log_q_ratio += (q_cur_given_prime - q_prime_given_cur)

    #### Update intermediate arrays ####

      data_cur = deepcopy(data_prime)
      track_cur = deepcopy(track_prime)

  end # end of for Δi in 1:size(Δ_array, 1)


  AddRem_track = [-1, 0, 0, -1, -1, -1, -1, -1]

  return(data_prime, track_prime, init_cond_cur, log_q_ratio, AddRem_track)
end



####################
### MH FUNCTIONS ###
####################

### Metropolis Hastings functions ####

function mult_MH_accept(post, post_prime, log_q_ratio, log_params_cur, log_params_draw)
  alpha =  (sum(log_params_draw) + post_prime) - (sum(log_params_cur) + post) + log_q_ratio

  return( min(0, alpha) )
end

function mh_accept_ratio(post, post_prime, log_q_ratio)

  log_α_ratio =  min( 0,  (post_prime) - (post) + log_q_ratio )

  return(log_α_ratio)
end

# MH function for parameters
function metropolis_hastings_step_params(data_cur, track_cur,
                                          post_cur, prior_dists,
                                          N_its, res, other_res,
                                          it, log_params_cur,
                                          n_tune, m, λ, init_cond,
                                          covarM, timestep_size)

  # Propose an update
  log_params_draw, log_q_ratio, mixture, λ_new = propose_params(N_its, res, other_res,
                                                                it, log_params_cur,
                                                                n_tune, m, λ, covarM)

  # Early return: Update is invalid
  if isfinite(log_q_ratio) == false
                    # is_accepted, log_α_ratio, post_cur, post_prime
    return(log_params_cur, [false, -Inf, post_cur, -Inf], mixture, λ_new)
  end

  # Calculate new llh arrays and posteriors
  post_prime = posterior(data_cur, track_cur, exp.(log_params_draw), prior_dists, init_cond, timestep_size)


  # Calculate MH acceptance probability
  log_α_ratio = mult_MH_accept(post_cur, post_prime, log_q_ratio, log_params_cur, log_params_draw)
  is_accepted = false

  # Accept/Reject
  accept_test = log(rand())
  if accept_test < log_α_ratio  #If yes:

    post_cur = deepcopy(post_prime)
    log_params_cur = deepcopy(log_params_draw)

    is_accepted = true
  end

  return(log_params_cur, post_cur, [is_accepted, log_α_ratio, post_cur, post_prime], mixture, λ_new)
end

# MH function for general data augmentation
function metropolis_hastings_step_aug(proposal_func, data_cur, track_cur, post_cur, prior_dists, log_params_cur, event, init_cond_cur, T, nt_up, timestep_size)

  # Propose an update
  data_prime, track_prime, init_cond_prime, log_q_ratio, update_tracker = proposal_func(data_cur, track_cur, event, exp.(log_params_cur), T, init_cond_cur, nt_up, timestep_size)

  # Early return: Update is invalid
  if isfinite(log_q_ratio) == false
                                  # is_accepted, log_α_ratio, post_cur, post_prime
    return(data_cur, track_cur, init_cond_cur, post_cur, [false, -Inf, post_cur, -Inf], update_tracker)
  end

  # Calculate new posteriors
  post_prime = posterior(data_prime, track_prime, exp.(log_params_cur), prior_dists, init_cond_prime, timestep_size)

  # Calculate MH acceptance probability
  log_α_ratio = mh_accept_ratio(post_cur, post_prime, log_q_ratio)
  is_accepted = false

  # Accept/Reject
  accept_test = log(rand())
  if accept_test < log_α_ratio  #If yes:

    data_cur = deepcopy(data_prime)
    track_cur = deepcopy(track_prime)
    init_cond_cur = deepcopy(init_cond_prime)

    post_cur = deepcopy(post_prime)

    update_tracker[2:3] = [1.,1.]

    is_accepted = true
  end

  return(data_cur, track_cur, init_cond_cur, post_cur, [is_accepted, log_α_ratio, post_cur, post_prime], update_tracker)
end


#####################
### Adaptive MCMC ###
#####################

function Initialise(data_in, track_in, params_true, prior_dists, init_conditions, timestep_size)

  # Current values
  β_cur, γ_cur, δ_cur = params_true
  params_draw = 0

  # True values
  β_true, γ_true, δ_true = params_true

  # Prior distributions
  d_β, d_γ, d_δ = prior_dists


  arb = true
  while arb == true
    # Initialise parameters
    if ismissing(β_true)
      β_cur = rand(d_β)
    end
    if ismissing(γ_true)
      γ_cur = rand(d_γ)
    end
    if ismissing(δ_true)
      δ_cur = rand(d_δ)
    end

    # Current draws
    params_draw = [β_cur, γ_cur, δ_cur]

    #Calculate the log-likelihood

    init_llh = llh(data_in, track_in, params_draw, init_conditions, timestep_size)

    if init_llh != -Inf
      break
    end
  end

  println("<", "Initialised!", ">")

  return(params_draw)
end

function create_results_arrays(N_its, its_per_frame)

  res = Array{Float64}(undef, N_its, 4)
  # :β, :γ, :δ :sample

  other_res = Array{Float64}(undef, N_its, 7)
  # :λ, :m, :acc, :log_α_ratio,
  # :post, :post_prime, :sample

  aug_res = Array{Float64}(undef, N_its, 16)
  # [is_accepted, log_α_ratio, post, post_prime] for
  # Move SE; Move EI;
  # Add/Rem SE; Add/Rem EI;

  move_SE_tracker = Array{Float64}(undef, N_its, 5)
  # :t, :is_accepted, :reason,
  # :Δ_time, :num_moved

  move_EI_tracker = Array{Float64}(undef, N_its, 5)
  # :t, :is_accepted, :reason,
  # :Δ_time, :num_moved

  AddRem_SE_tracker = Array{Float64}(undef, N_its, 8)
  # :t, :is_accepted, :reason,
  # :SE_before, :SE_after, :Δ_diff, :cS, :prob

  AddRem_EI_tracker = Array{Float64}(undef, N_its, 8)
  # :t, :is_accepted, :reason,
  # :EI_before, :EI_after, :Δ_diff, :cE, :prob

  no_of_saves = convert(Int64, N_its/its_per_frame)
  SIM_tracker = Array{Array{}}(undef, no_of_saves)
  EVENTS_tracker = Array{Array{}}(undef, no_of_saves)

  return(res, other_res, aug_res, move_SE_tracker, move_EI_tracker, AddRem_SE_tracker, AddRem_EI_tracker, SIM_tracker, EVENTS_tracker)
end

function rename_results_arrays(res, other_res, aug_res, move_SE_tracker, move_EI_tracker, AddRem_SE_tracker, AddRem_EI_tracker, SIM_tracker, EVENTS_tracker)

  res[:, 1:3] = exp.(res[:, 1:3])

  res = DataFrame(res, :auto)
  rename!(res, [:β, :γ, :δ, :sample])

  other_res = DataFrame(other_res, :auto)
  rename!(other_res, [:λ, :m, :acc, :log_α_ratio,
                      :post, :post_prime, :sample])

  aug_res = DataFrame(aug_res, :auto)
  rename!(aug_res, [:is_accepted_move_SE, :log_α_ratio_move_SE, :post_move_SE, :post_prime_move_SE,
                      :is_accepted_move_EI, :log_α_ratio_move_EI, :post_move_EI, :post_prime_move_EI,
                      :is_accepted_AddRem_SE, :log_α_ratio_AddRem_SE, :post_AddRem_SE, :post_prime_AddRem_SE,
                      :is_accepted_AddRem_EI, :log_α_ratio_AddRem_EI, :post_AddRem_EI, :post_prime_AddRem_EI])

  move_SE_tracker = DataFrame(move_SE_tracker, :auto)
  rename!(move_SE_tracker, [:t, :is_accepted, :reason,
                            :Δ_time, :num_moved])

  move_EI_tracker = DataFrame(move_EI_tracker, :auto)
  rename!(move_EI_tracker, [:t, :is_accepted, :reason,
                            :Δ_time, :num_moved])

  AddRem_SE_tracker = DataFrame(AddRem_SE_tracker, :auto)
  rename!(AddRem_SE_tracker, [:t, :is_accepted, :reason,
                              :SE_before, :SE_after, :Δ_diff, :cS, :prob])

  AddRem_EI_tracker = DataFrame(AddRem_EI_tracker, :auto)
  rename!(AddRem_EI_tracker, [:t, :is_accepted, :reason,
                              :EI_before, :EI_after, :Δ_diff, :cE, :prob])

  return(res, other_res, aug_res, move_SE_tracker, move_EI_tracker, AddRem_SE_tracker, AddRem_EI_tracker, SIM_tracker, EVENTS_tracker)
end

function Blk_Adaptive_MCMC_post_all(;N_its, data_init, track_init,
                                    prior_dists, params_true,
                                    λ_init, m_init,
                                    data_aug_infer, param_infer,
                                    init_cond, T, nt_up, timestep_size)

  ##############################################
  ### Initialise the parameters and epidemic ###
  ##############################################

      params_cur = Initialise(data_init, track_init, params_true, prior_dists, init_cond, timestep_size)

      log_params_cur = log.(params_cur)

      covarM = CovMatrix()

  ########################
  ### Results matrices ###
  ########################

    res, other_res, aug_res,
     move_SE_tracker, move_EI_tracker,
      AddRem_SE_tracker, AddRem_EI_tracker,
       SIM_tracker, EVENTS_tracker = create_results_arrays(N_its, 100)

  ##########################
  ### Functional objects ###
  ##########################

      it = 1

      λ = λ_init
      m = m_init
      Δm = m_init/100

      init_cond_cur = init_cond
      data_cur = data_init
      track_cur = track_init

      post_cur = posterior(data_cur, track_cur, params_cur, prior_dists, init_cond_cur, timestep_size)

      n_tune = 1

  ###########################
  ### ~~ THE ALGORITHM ~~ ###
  ###########################

      while it <= N_its

        ##############################
        ### MH step for parameters ###
        ##############################

        mh_res = [-Inf, -Inf, Inf, Inf]

        if param_infer == true

          if it > 1
            fit!(covarM, res[(it-1), 1:3])
          end

          log_params_cur, post_cur, mh_res, mixture, λ =
                    metropolis_hastings_step_params(data_cur, track_cur,
                                                    post_cur, prior_dists,
                                                    N_its, res, other_res,
                                                    it, log_params_cur,
                                                    n_tune, m, λ, init_cond_cur,
                                                    covarM, timestep_size)


          if mixture >= 0.05 # mixture is set to 0 while tuning λ
            if mh_res[1] == 0
              m = m - (Δm/((it)^0.5))
            else
              m = m + 2.3*(Δm/((it)^0.5))
            end
          end

        end


        ###############################
        ### Data Augmentation Steps ###
        ###############################

        ###########################
        ### Move S→E Init Conds ###
        ###########################

        mh_res_init_conds_SE = [-Inf, -Inf, Inf, Inf]
        init_conds_SE_track = [-Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[1] == true

          data_cur, track_cur, init_cond_cur, post_cur, mh_res_init_conds_SE, init_conds_SE_track =
                metropolis_hastings_step_aug(propose_init, data_cur, track_cur, post_cur, prior_dists, log_params_cur, 1, init_cond_cur, T, nt_up[1], timestep_size)

        end


        ###########################
        ### Move E→I Init Conds ###
        ###########################

        mh_res_init_conds_EI = [-Inf, -Inf, Inf, Inf]
        init_conds_EI_track = [-Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[2] == true

          data_cur, track_cur, init_cond_cur, post_cur, mh_res_init_conds_EI, init_conds_EI_track =
                metropolis_hastings_step_aug(propose_init, data_cur, track_cur, post_cur, prior_dists, log_params_cur, 2, init_cond_cur, T, nt_up[1], timestep_size)

        end


        ###################################
        ### Move S→E event through time ###
        ###################################

        mh_res_move_SE = [-Inf, -Inf, Inf, Inf]
        move_SE_track = [-Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[3] == true

          data_cur, track_cur, init_cond_cur, post_cur, mh_res_move_SE, move_SE_track =
                metropolis_hastings_step_aug(propose_move, data_cur, track_cur, post_cur, prior_dists, log_params_cur, 1, init_cond_cur, T, nt_up[1], timestep_size)

        end

        ###################################
        ### Move E→I event through time ###
        ###################################

        mh_res_move_EI = [-Inf, -Inf, Inf, Inf]
        move_EI_track = [-Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[4] == true

          data_cur, track_cur, init_cond_cur, post_cur, mh_res_move_EI, move_EI_track =
                metropolis_hastings_step_aug(propose_move, data_cur, track_cur, post_cur, prior_dists, log_params_cur, 2, init_cond_cur, T, nt_up[1], timestep_size)


        end

        ###################################
        ### AddRem SE event through time ##
        ###################################

        mh_res_AddRem_SE = [-Inf, -Inf, Inf, Inf]
        AddRem_SE_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[5] == true

          data_cur, track_cur, init_cond_cur, post_cur, mh_res_AddRem_SE, AddRem_SE_track =
                metropolis_hastings_step_aug(propose_AddRem, data_cur, track_cur, post_cur, prior_dists, log_params_cur, 1, init_cond_cur, T, nt_up[2], timestep_size)

        end

        ###################################
        ### AddRem EI event through time ##
        ###################################

        mh_res_AddRem_EI = [-Inf, -Inf, Inf, Inf]
        AddRem_EI_track = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]

        if data_aug_infer[6] == true

          data_cur, track_cur, init_cond_cur, post_cur, mh_res_AddRem_EI, AddRem_EI_track =
                metropolis_hastings_step_aug(propose_AddRem, data_cur, track_cur, post_cur, prior_dists, log_params_cur, 2, init_cond_cur, T, nt_up[2], timestep_size)

        end

        ##########################
        ### Record the results ###
        ##########################

          # Record parameters
          res[it,:] = [log_params_cur; it]

          # Record other
          other_res[it,:] = [λ, m, mh_res[1], mh_res[2], mh_res[3], mh_res[4], it]

          # Record data aug
          aug_res[it, :] = [mh_res_move_SE ; mh_res_move_EI ; mh_res_AddRem_SE ; mh_res_AddRem_EI]

          # Record update tracking data
          move_SE_tracker[it, :] = move_SE_track
          move_EI_tracker[it, :] = move_EI_track
          AddRem_SE_tracker[it, :] = AddRem_SE_track
          AddRem_EI_tracker[it, :] = AddRem_EI_track

          if rem(it, 100) == 0
            SIM_count = convert(Int64, it/100)
            SIM_tracker[SIM_count] = data_cur
            EVENTS_tracker[SIM_count] = track_cur
            println("it = ", it)
          end

          # Update count
          if it == (n_tune*25)
            n_tune = n_tune + 1
          end
          it = it + 1

      end #end of while

  res, other_res, aug_res, move_SE_tracker, move_EI_tracker, AddRem_SE_tracker, AddRem_EI_tracker, SIM_tracker, EVENTS_tracker = rename_results_arrays(res, other_res, aug_res, move_SE_tracker, move_EI_tracker, AddRem_SE_tracker, AddRem_EI_tracker, SIM_tracker, EVENTS_tracker)

  return(res, other_res, aug_res, move_SE_tracker, move_EI_tracker, AddRem_SE_tracker, AddRem_EI_tracker, SIM_tracker, EVENTS_tracker)
end
