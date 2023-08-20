# Required packages
using Random
using Distributions
using DataFrames
using Plots
using CSV

# Function to perform a single Gillespie simulation step
function gillespie_step(S::Int, E::Int, I::Int, R::Int, β::Float64, δ::Float64, γ::Float64)
    # Calculate total population and rates
    N = S + E + I + R
    λ_t = β * S * I/N
    δ_t = δ * E
    γ_t = γ * I

    events = [0,0,0]

    # Calculate total rate of events
    total_rate = λ_t + δ_t + γ_t

    # Time until the next event occurs
    τ = rand(Exponential(1/total_rate))

    # Randomly choose which event occurs
    rand_val = rand()
    if rand_val < λ_t / total_rate
        S -= 1
        E += 1
        events += [1,0,0]
    elseif rand_val < (λ_t + δ_t) / total_rate
        E -= 1
        I += 1
        events += [0,1,0]
    else
        I -= 1
        R += 1
        events += [0,0,1]
    end

    return S, E, I, R, τ, events
end

function run_seir_simulation(β::Float64, δ::Float64, γ::Float64, initial_S::Int, initial_E::Int, initial_I::Int, initial_R::Int, max_steps::Int = 10_000)
    # Initialize the compartments and time
    S = initial_S
    E = initial_E
    I = initial_I
    R = initial_R
    t = 0.0

    # Pre-allocate a large data frame to store the simulation results
    df = DataFrame(
        Time = zeros(Float64, max_steps),
        Susceptible = zeros(Int, max_steps),
        Exposed = zeros(Int, max_steps),
        Infected = zeros(Int, max_steps),
        Recovered = zeros(Int, max_steps),
        SE_event = zeros(Int, max_steps),
        EI_event = zeros(Int, max_steps),
        IR_event = zeros(Int, max_steps)
    )

    # Count variable to track the current index
    count = 1
    events = [0,0,0]

    # Run the simulation until there are no more exposed or infectious individuals
    while E > 0 || I > 0
        # Add the current values to the data frame
        df.Time[count] = t
        df.Susceptible[count] = S
        df.Exposed[count] = E
        df.Infected[count] = I
        df.Recovered[count] = R
        df[count, 6] += events[1]
        df[count, 7] += events[2]
        df[count, 8] += events[3]

        S, E, I, R, τ, events = gillespie_step(S, E, I, R, β, δ, γ)

        # Update time
        t += τ

        # Increment the count variable
        count += 1

        # Check if the data frame is full
        if count > max_steps
            throw(ArgumentError("Simulation ran for too many steps, increase max_steps parameter"))
        end
    end

    # Add the current values to the data frame
    df.Time[count] = t
    df.Susceptible[count] = S
    df.Exposed[count] = E
    df.Infected[count] = I
    df.Recovered[count] = R
    df[count, 6] += events[1]
    df[count, 7] += events[2]
    df[count, 8] += events[3]

    # Remove unused rows from the data frame
    df = df[1:count, :]

    return df
end

function aggregate_simulation_data(df::DataFrame, time_step::Float64)

    # Create an empty data frame to store the aggregated results
    aggregated_df = DataFrame(
        Time = Float64[],
        Susceptible = Int[],
        Exposed = Int[],
        Infected = Int[],
        Recovered = Int[],
        SE_events = Int[],
        EI_events = Int[],
        IR_events = Int[]
    )

    # Create timestep array
    max_time = maximum(df[:, 1])
    max_timestep = Int(ceil((ceil(max_time/time_step)) * time_step))

    seq = collect(0:time_step:max_timestep)


    # Initialize variables to track population counts
    S_count = df.Susceptible[1]
    E_count = df.Exposed[1]
    I_count = df.Infected[1]
    R_count = df.Recovered[1]

    # Fill in states and events
    for i in 2:length(seq)

        println(seq[i])

        rows = seq[i-1] .< df.Time .<= seq[i]

        df_timestep = df[rows, :]

        SE_events = sum(df_timestep[:, 6])
        EI_events = sum(df_timestep[:, 7])
        IR_events = sum(df_timestep[:, 8])

        S_count = S_count - SE_events
        E_count = E_count + SE_events - EI_events
        I_count = I_count + EI_events - IR_events
        R_count = R_count + IR_events


        push!(aggregated_df, (seq[i], S_count, E_count, I_count, R_count, SE_events, EI_events, IR_events))

    end

    return(aggregated_df)
end

function generate_EI_and_SE_events(IR_t, β, δ, Δt)

    p_exp = 1 - exp(-β)

    p_inf = 1-exp(-δ)

    EI_t = -10000
    SE_t = -10000

    while EI_t == SE_t

        EI_t = IR_t - 1*Δt - rand(Geometric(p_inf))*Δt
        SE_t = EI_t - 1*Δt - rand(Geometric(p_exp))*Δt

    end

    SE_t = round(SE_t, digits = 5)
    EI_t = round(EI_t, digits = 5)

    return(SE_t, EI_t)
end

function generate_valid_epidemic(agg_df, β, δ, Δt)

    # Create an empty data frame to store the aggregated results
    valid_agg_df = DataFrame(
        Time = agg_df.Time,
        Susceptible = zeros(Int, length(agg_df.Time)),
        Exposed = zeros(Int, length(agg_df.Time)),
        Infected = zeros(Int, length(agg_df.Time)),
        Recovered = agg_df.Recovered,
        SE_events = zeros(Int, length(agg_df.Time)),
        EI_events = zeros(Int, length(agg_df.Time)),
        IR_events = agg_df.IR_events
    )

    seq = Array(agg_df.Time)

    # Set initial conditions
    new_init_conds = [sum(agg_df[1, 2:5]), 0, 0, 0]
    first = 1


    # Generate new SE and EI events
    for (t, num_R) in zip(agg_df.Time, agg_df.IR_events)
        println(t, "   ", num_R)

        if first == 1
            if num_R > 0
                new_init_conds[1] -= 1
                new_init_conds[3] += 1
                num_R -= 1
                first = 0
            end
        end

        if num_R > 0
            for i in 1:num_R

                # Generate new event times
                SE_t, EI_t = generate_EI_and_SE_events(t, β, δ, Δt)
                println([SE_t, EI_t])

                # If before 0, add to initial conditions
                if SE_t <= 0

                    new_init_conds[1] -= 1
                    new_init_conds[2] += 1

                end

                if EI_t <= 0

                    new_init_conds[2] -= 1
                    new_init_conds[3] += 1

                end


                # If during epidemic, apply to appropriate timestep
                if SE_t > 0

                    pos = findfirst(isequal(SE_t), seq)

                    valid_agg_df[pos:end, 2] .-= 1
                    valid_agg_df[pos:end, 3] .+= 1

                    valid_agg_df[pos, 6] += 1

                end

                if EI_t > 0

                    pos = findfirst(isequal(EI_t), seq)

                    valid_agg_df[pos:end, 3] .-= 1
                    valid_agg_df[pos:end, 4] .+= 1

                    valid_agg_df[pos, 7] += 1

                end

                pos = findfirst(isequal(t), seq)

                valid_agg_df[pos:end, 4] .-= 1

            end
        end
    end

    valid_agg_df[:, 2] .+= new_init_conds[1]
    valid_agg_df[:, 3] .+= new_init_conds[2]
    valid_agg_df[:, 4] .+= new_init_conds[3]
    valid_agg_df[:, 5] .+= new_init_conds[4]

    return(new_init_conds, valid_agg_df)
end



# SEIR parameters
β = 0.25     # Infection rate
δ = 0.08
γ = 0.22   # Recovery rate

# Example usage
initial_S = 999
initial_E = 0
initial_I = 1
initial_R = 0
max_steps = 10000

Random.seed!(38)

simulation_df = run_seir_simulation(β, δ, γ, initial_S, initial_E, initial_I, initial_R, max_steps)

CSV.write("Data/simulated_data_1000_continuous.csv", simulation_df, header = false)

# println(simulation_df[end, :])

# Plotting
# plot(simulation_df.Time, simulation_df.Susceptible, label="Susceptible", color=:green, seriestype=:steppre, xlabel="Time", ylabel="Population Size", linewidth=1)
#     plot!(simulation_df.Time, simulation_df.Exposed, label="Exposed", color=:yellow, seriestype=:steppre, linewidth=1)
#     plot!(simulation_df.Time, simulation_df.Infected, label="Infected", color=:red, seriestype=:steppre, linewidth=1)
#     plot!(simulation_df.Time, simulation_df.Recovered, label="Recovered", color=:blue, seriestype=:steppre, linewidth=1)
#     plot!(legend=:outertopright)

ongoing_sim = simulation_df[simulation_df[:,5] .<= 250, :]

CSV.write("Data/simulated_data_1000_continuous_cut.csv", ongoing_sim, header = false)


sort(copy((ongoing_sim[2:end, 1] - ongoing_sim[1:(end-1), 1])))[1:15]

plot(ongoing_sim.Time, ongoing_sim.Susceptible, label="Susceptible", color=:green, seriestype=:steppre, xlabel="Time", ylabel="Population Size", linewidth=1)
    plot!(ongoing_sim.Time, ongoing_sim.Exposed, label="Exposed", color=:yellow, seriestype=:steppre, linewidth=1)
    plot!(ongoing_sim.Time, ongoing_sim.Infected, label="Infected", color=:red, seriestype=:steppre, linewidth=1)
    plot!(ongoing_sim.Time, ongoing_sim.Recovered, label="Recovered", color=:blue, seriestype=:steppre, linewidth=1)
    plot!(legend=:outertopright)



### AGGREGATION ####

aggregated_df_point2 = aggregate_simulation_data(ongoing_sim, 0.2)

valid_init_conds_point2, valid_agg_sim_point2 = generate_valid_epidemic(aggregated_df_point2, 0.25, 0.08, 0.2)

# Plotting
plot(aggregated_df_point2.Time, aggregated_df_point2.Susceptible, label="Susceptible", color=:green, seriestype=:steppost, xlabel="Time", ylabel="Population Size", linewidth=1)
    plot!(aggregated_df_point2.Time, aggregated_df_point2.Exposed, label="Exposed", color=:yellow, seriestype=:steppost, linewidth=1)
    plot!(aggregated_df_point2.Time, aggregated_df_point2.Infected, label="Infected", color=:red, seriestype=:steppost, linewidth=1)
    plot!(aggregated_df_point2.Time, aggregated_df_point2.Recovered, label="Recovered", color=:blue, seriestype=:steppost, linewidth=1, ylims = [0, 1000])
    plot!(legend=:outertopright)


CSV.write("Data/init_conds_agg_point2.csv", Tables.table(valid_init_conds_point2), header = false)
CSV.write("Data/simulated_data_1000_agg_point2.csv", valid_agg_sim_point2[:, 2:5], header = false)
CSV.write("Data/simulated_track_1000_agg_point2.csv", valid_agg_sim_point2[:, 6:8], header = false)

full_aggregated_df_point2 = aggregate_simulation_data(simulation_df, 0.2)
CSV.write("Data/full_agg_point2.csv", full_aggregated_df_point2[:, 2:5], header = false)


### AGGREGATION ####

aggregated_df_1 = aggregate_simulation_data(ongoing_sim, 1.0)

valid_init_conds_1, valid_agg_sim_1 = generate_valid_epidemic(aggregated_df_1, 0.25, 0.08, 1)


# Plotting
plot(aggregated_df_1.Time, aggregated_df_1.Susceptible, label="Susceptible", color=:green, seriestype=:steppost, xlabel="Time", ylabel="Population Size", linewidth=1)
    plot!(aggregated_df_1.Time, aggregated_df_1.Exposed, label="Exposed", color=:yellow, seriestype=:steppost, linewidth=1)
    plot!(aggregated_df_1.Time, aggregated_df_1.Infected, label="Infected", color=:red, seriestype=:steppost, linewidth=1)
    plot!(aggregated_df_1.Time, aggregated_df_1.Recovered, label="Recovered", color=:blue, seriestype=:steppost, linewidth=1, ylims = [0, 1000])
    plot!(legend=:outertopright)


CSV.write("Data/init_conds_agg_1.csv", Tables.table(valid_init_conds_1), header = false)
CSV.write("Data/simulated_data_1000_agg_1.csv", valid_agg_sim_1[:, 2:5], header = false)
CSV.write("Data/simulated_track_1000_agg_1.csv", valid_agg_sim_1[:, 6:8], header = false)

full_aggregated_df_1 = aggregate_simulation_data(simulation_df, 1.0)
CSV.write("Data/full_agg_1.csv", full_aggregated_df_1[:, 2:5], header = false)


### AGGREGATION ####

aggregated_df_7 = aggregate_simulation_data(ongoing_sim, 7.0)

valid_init_conds_7, valid_agg_sim_7 = generate_valid_epidemic(aggregated_df_7, 0.25, 0.08, 7)


# Plotting
plot(aggregated_df_7.Time, aggregated_df_7.Susceptible, label="Susceptible", color=:green, seriestype=:steppost, xlabel="Time", ylabel="Population Size", linewidth=1)
    plot!(aggregated_df_7.Time, aggregated_df_7.Exposed, label="Exposed", color=:yellow, seriestype=:steppost, linewidth=1)
    plot!(aggregated_df_7.Time, aggregated_df_7.Infected, label="Infected", color=:red, seriestype=:steppost, linewidth=1)
    plot!(aggregated_df_7.Time, aggregated_df_7.Recovered, label="Recovered", color=:blue, seriestype=:steppost, linewidth=1, ylims = [0, 1000])
    plot!(legend=:outertopright)


CSV.write("Data/init_conds_agg_7.csv", Tables.table(valid_init_conds_7), header = false)
CSV.write("Data/simulated_data_1000_agg_7.csv", valid_agg_sim_7[:, 2:5], header = false)
CSV.write("Data/simulated_track_1000_agg_7.csv", valid_agg_sim_7[:, 6:8], header = false)


full_aggregated_df_7 = aggregate_simulation_data(simulation_df, 7.0)
CSV.write("Data/full_agg_7.csv", full_aggregated_df_7[:, 2:5], header = false)



### AGGREGATION ####

aggregated_df_30 = aggregate_simulation_data(ongoing_sim, 30.0)

valid_init_conds_30, valid_agg_sim_30 = generate_valid_epidemic(aggregated_df_30, 0.25, 0.08, 30)


# Plotting
plot(aggregated_df_30.Time, aggregated_df_30.Susceptible, label="Susceptible", color=:green, seriestype=:steppost, xlabel="Time", ylabel="Population Size", linewidth=1)
    plot!(aggregated_df_30.Time, aggregated_df_30.Exposed, label="Exposed", color=:yellow, seriestype=:steppost, linewidth=1)
    plot!(aggregated_df_30.Time, aggregated_df_30.Infected, label="Infected", color=:red, seriestype=:steppost, linewidth=1)
    plot!(aggregated_df_30.Time, aggregated_df_30.Recovered, label="Recovered", color=:blue, seriestype=:steppost, linewidth=1, ylims = [0, 1000])
    plot!(legend=:outertopright)


CSV.write("Data/init_conds_agg_30.csv", Tables.table(valid_init_conds_30), header = false)
CSV.write("Data/simulated_data_1000_agg_30.csv", valid_agg_sim_30[:, 2:5], header = false)
CSV.write("Data/simulated_track_1000_agg_30.csv", valid_agg_sim_30[:, 6:8], header = false)

full_aggregated_df_30 = aggregate_simulation_data(simulation_df, 30.0)
CSV.write("Data/full_agg_30.csv", full_aggregated_df_30[:, 2:5], header = false)
