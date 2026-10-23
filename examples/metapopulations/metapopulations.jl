# Run from base jOpqua directory as
# julia --project=. examples/metapopulations/metapopulations.jl
# unless viewing flamegraph, then run from console

# using Revise
using jOpqua

using StaticArrays
using Random

using Distances

using BenchmarkTools
using ProfileView

# Model setup
function run(seed::Int64, t_vec::Vector{Float64})
    # Parameters

    start_genome = "AAAA"

    pat_type = jOpqua.newPathogenType(
        "pat_type",
        num_loci=4,
        possible_alleles="AB",
        receiveContactHostwideCoefficient=(s::String,pop_id::String) -> 0.0,
    )

    hos_type = jOpqua.newHostType("hos_type")

    pop_type = jOpqua.newPopulationType(
        "pop_type",
        clearance_coefficient=1.0,
        contact_coefficient=1.1,
        receive_contact_coefficient=1.0,
        inoculum_coefficient=1.0,
        transition_coefficient=1.0e-3,
        # pathogenFractions=jOpqua.pathogenFractionsProportionalFitness,
    )

    num_hosts = 1000
    num_infected = Int(num_hosts * 0.1)
    host_genome = ""

    # Setup
    model = jOpqua.newModel()
    pop1 = jOpqua.newPopulation!("pop1", pop_type, model)
    pop2 = jOpqua.newPopulation!("pop2", pop_type, model)
    jOpqua.setPopulationContactCoefficient!(pop1, pop2, 1.0, model)
    jOpqua.setPopulationTransitionCoefficient!(pop2, pop1, 1.0, model)
    jOpqua.addHostsToPopulation!(num_hosts, host_genome, hos_type, pop1, model)
    jOpqua.addHostsToPopulation!(num_hosts, host_genome, hos_type, pop2, model)
    pat = jOpqua.newPathogen!(start_genome, pop1, pat_type)

    for h in 1:num_infected
        jOpqua.addPathogenToHost!(pat, h, pop1, model)
    end

    # Simulate
    Random.seed!(seed)
    model, output = jOpqua.simulate!(
        model, t_vec, population_host_samples=Dict("pop1" => 200, "pop2" => 200)
    )

    # Data output and plots
    compartment_data = jOpqua.saveCompartments(output, "examples/metapopulations/compartment_pathogen_evolution.csv")
    jOpqua.plotCompartments(compartment_data, ["pop1"], "examples/metapopulations/compartment_pathogen_evolution_pop1.png")
    jOpqua.plotCompartments(compartment_data, ["pop2"], "examples/metapopulations/compartment_pathogen_evolution_pop2.png")
end

run(1, collect(0.0:2.0:4.0)) # compile
@time run(0, collect(0.0:2.0:1500.0))
# Total events: 1728203
#   4.460974 seconds (55.20 M allocations: 7.560 GiB, 15.71% gc time, 1.58% compilation time)
# 25 May 2026 Julia 1.12.6 Apple M3 Max 128 GB RAM

# @profview run(2, collect(0.0:2.0:1500.0))
