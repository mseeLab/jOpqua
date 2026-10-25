# Run from base jOpqua directory as
# julia --project=. examples/acquired_immunity/acquired_immunity.jl
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
    optimal_genome = "BBBB"
    genome_length = length(optimal_genome)

    pat_type = jOpqua.newPathogenType(
        "pat_type",
        num_loci=genome_length,
        possible_alleles="AB",
        contactSpecificCoefficient=(s::String,pop_id::String) -> 1.0 + (0.1 * (genome_length - hamming(s, optimal_genome)) / genome_length),
        receiveContactHostwideCoefficient=(s::String,pop_id::String) -> 0.0, # makes infected hosts immune to superinfection
    )

    res_type = jOpqua.newResponseType(
        "res_type",
        reactivityCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String, pop_id::String) -> imp_g == pat_g ? 1.0 : 0.0,
        # transmissionEfficiencyInteractionSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> imp_g == pat_g ? 0.0 : 1.0,
        clearanceInteractionSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String, pop_id::String) -> imp_g == pat_g ? 1.1e0 : 1.0,
        # responseAcquisitionInteractionSpecificCoefficient=(hos_g::String, imp_g::String, mat_g::String, pat_g::String) -> imp_g == pat_g ? 0.0 : 1.0,
    )

    hos_type = jOpqua.newHostType("hos_type")

    pop_type = jOpqua.newPopulationType(
        "pop_type",
        clearance_coefficient=1.0,
        contact_coefficient=1.05,
        # response_acquisition_coefficient=0.5,
        response_acquisition_upon_clearance_coefficient=0.1,
        response_loss_coefficient=0.01,
        receive_contact_coefficient=1.0,
        mutations_per_generation_coefficient=0.0005,
        inoculum_coefficient=1.0,
        # death_coefficient=0.001,
        # birth_coefficient=0.001,
        pathogenFractions=jOpqua.pathogenFractionsProportionalFitness,
        response_types=Dict{String,jOpqua.ResponseType}([(res_type.id => res_type)]),
        developResponses=(
            pathogen::jOpqua.Pathogen, host::jOpqua.Host, pop_id::String,
            existing_responses::Dict{Tuple{String,String,String,String},jOpqua.Response},
            response_types::Dict{String,jOpqua.ResponseType},
            birth_time::Float64
        ) -> [jOpqua.deNovoResponse(
            pathogen, host, "pop_type", existing_responses, response_types, birth_time;
            response_type_id="res_type"
        )],
    )

    num_hosts = 10000
    num_infected = Int(num_hosts * 0.05)
    host_genome = ""

    # Setup
    model = jOpqua.newModel()
    pop = jOpqua.newPopulation!("pop", pop_type, model)
    jOpqua.addHostsToPopulation!(num_hosts, host_genome, hos_type, pop, model)
    pat = jOpqua.newPathogen!(start_genome, pop, pat_type)

    for h in 1:num_infected
        jOpqua.addPathogenToHost!(pat, h, pop, model)
    end

    # Simulate
    Random.seed!(seed)
    model, output = jOpqua.simulate!(
        model, t_vec, population_host_samples=Dict("pop" => 200)
    )

    # Data output and plots
    compartment_data = jOpqua.saveCompartments(output, "examples/acquired_immunity/compartment_acquired_immunity.csv")
    jOpqua.plotCompartments(compartment_data, ["pop"], "examples/acquired_immunity/compartment_acquired_immunity.png")

    his_dat = jOpqua.saveHistory(output, "examples/acquired_immunity/history_acquired_immunity.csv")
    composition_data = jOpqua.saveComposition(
        his_dat, "examples/acquired_immunity/composition_acquired_immunity.csv",
        num_top_sequences=7, #track_specific_sequences=["AAAA", "BBBB"]
    )
    jOpqua.plotComposition(
        composition_data, "examples/acquired_immunity/composition_acquired_immunity.png",
        normalized=true, ylabel="Fraction",
    )

    composition_data = jOpqua.saveComposition(
        his_dat, "examples/acquired_immunity/composition_acquired_immunity_responses.csv",
        num_top_sequences=7, #track_specific_sequences=["AAAA", "BBBB"],
        type_of_composition="Response_imprinted_sequence"
    )
    jOpqua.plotComposition(
        composition_data, "examples/acquired_immunity/composition_acquired_immunity_responses.png",
        # normalized=true, ylabel="Fraction",
        normalized=false, ylabel="Number",
    )

    # nwks = jOpqua.saveNewick(output, "examples/acquired_immunity/pathogen_newick_acquired_immunity.nwk", branch_length="Time", info_separator=" ")
    # for nwk in nwks
    #     jOpqua.plotPhylogeny(nwk, "examples/acquired_immunity/pathogen_newick_acquired_immunity.png")
    # end
end

run(1, collect(0.0:2.0:4.0)) # compile
@time run(10, collect(0.0:2.0:1500.0))
# Total events: 3342817
#   7.751170 seconds (82.72 M allocations: 14.894 GiB, 20.86% gc time, 7.99% compilation time: <1% of which was recompilation)
# 25 May 2026 Julia 1.12.6 Apple M3 Max 128 GB RAM

# @profview run(2, collect(0.0:2.0:1500.0))
