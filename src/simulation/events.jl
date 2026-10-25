using StaticArrays
using StatsBase
using PoissonRandom
using Random
using Flexle

# General actions

function mutantSequence!(
    sequence::String, num_loci::Int64, possible_alleles::String,
    mean_mutations_per_replication::Float64)

    chromosomes = split(sequence, CHROMOSOME_SEPARATOR)
    for (i, chromosome_pair) in enumerate(chromosomes)
        homologous_chromosomes = split(chromosome_pair, HOMOLOGOUS_CHROMOSOME_SEPARATOR)
        for (j, chromosome) in enumerate(homologous_chromosomes)
            x = zeroTruncatedPoisson(mean_mutations_per_replication * length(chromosome) / num_loci)
            loci = sort(unique(rand(
                1:length(chromosome), x
            )))
            seq = ""
            for locus in loci
                seq = seq * chromosome[length(seq)+1:locus-1] * rand(possible_alleles)
            end
            homologous_chromosomes[j] = seq * sequence[length(seq)+1:end]
        end
        chromosomes[i] = join(homologous_chromosomes, HOMOLOGOUS_CHROMOSOME_SEPARATOR)
    end

    return join(chromosomes, CHROMOSOME_SEPARATOR)
end

function mutantPathogen!(
        pathogen::Pathogen, host::Host, population::Population, birth_time::Float64, expected_mutations::Float64)
    seq = mutantSequence!(
        pathogen.sequence, pathogen.type.num_loci, pathogen.type.possible_alleles,
        expected_mutations
    )

    if haskey(population.pathogens, seq)
        return population.pathogens[seq]
    else
        return newPathogen!(
            seq, population, pathogen.type,
            parents=MVector{2,Union{Pathogen,Nothing}}([pathogen, nothing]),
            birth_time=birth_time
        )
    end
end

function recombinantSequences(
    seq_1::String, seq_2::String, num_loci::Int64,
    mean_recombination_crossovers_1::Float64, mean_recombination_crossovers_2::Float64)

    children = MVector{2,String}(seq_1, seq_2)

    if mean_recombination_crossovers_1 > 0.0 && mean_recombination_crossovers_2 > 0.0
        num_evts = zeroTruncatedPoisson((
            mean_recombination_crossovers_1 +
            mean_recombination_crossovers_2
        ) / 2.0 )
        loci = rand(1:num_loci, num_evts)
        for l in loci
            temp_child_1 = children[1]
            children = MVector{2,String}(children[1][1:l-1] * children[2][l:end], children[2])
            children = MVector{2,String}(children[1], children[2][1:l-1] * temp_child_1[l:end])
            # MVector of Strings could not be indexed-set
        end
    else
        num_evts = 0
    end

    children = MVector{2,Vector{String}}([split(seq, CHROMOSOME_SEPARATOR) for seq in children])
    parent = rand(0:1, length(children[1]))

    children = SVector{2,String}([
        join([children[parent[i]+1][i] for i in 1:length(children[1])], CHROMOSOME_SEPARATOR),
        join([children[(parent[i]!=true)+1][i] for i in 1:length(children[2])], CHROMOSOME_SEPARATOR)
    ])

    return children
end

function recombinantPathogens!(
        pathogen_1::Pathogen, pathogen_2::Pathogen, host::Host, population::Population,
        birth_time::Float64, mean_recombination_crossovers_1::Float64, mean_recombination_crossovers_2::Float64)
    children = recombinantSequences(
        pathogen_1.sequence, pathogen_2.sequence, pathogen_1.type.num_loci,
        mean_recombination_crossovers_1, mean_recombination_crossovers_2
    )

    if !haskey(population.pathogens, children[1])
        newPathogen!(
            children[1], population, pathogen_1.type,
            parents=MVector{2,Union{Pathogen,Nothing}}([pathogen_1, pathogen_2]),
            birth_time=birth_time
        )
    end

    return population.pathogens[children[1]]
end

function generateGamete(
    sequence::String, num_loci::Int64,
    possible_alleles::String, mean_mutations_per_replication::Float64)
    gamete = join(
        [
            rand(split(homologous_chromosome_pair, HOMOLOGOUS_CHROMOSOME_SEPARATOR))
            for homologous_chromosome_pair in split(sequence, CHROMOSOME_SEPARATOR)
        ],
        CHROMOSOME_SEPARATOR
    )
    if mean_mutations_per_replication > 0.0
        gamete = mutantSequence!(
            gamete, num_loci, possible_alleles,
            mean_mutations_per_replication
        )
    end

    return gamete
end

function generateZygote(gamete_1::String, gamete_2::String)
    chromosomes_1 = split(gamete_1, CHROMOSOME_SEPARATOR)
    chromosomes_2 = split(gamete_2, CHROMOSOME_SEPARATOR)
    zygote = ""
    for chromosome_idx in eachindex(chromosomes_1)
        zygote = zygote * chromosomes_1[chromosome_idx] * HOMOLOGOUS_CHROMOSOME_SEPARATOR * chromosomes_2[chromosome_idx] * CHROMOSOME_SEPARATOR
    end

    return zygote[1:end-length(CHROMOSOME_SEPARATOR)]
end

function addPathogenToHost!(pathogen::Pathogen, host_idx::Int64, population::Population, model::Model)
    if length(population.hosts[host_idx].pathogens) == 0
        if length(population.hosts[host_idx].responses) == 0
            population.compartment_vars[INFECTED_NAIVE] += 1
            population.compartment_vars[UNINFECTED_NAIVE] -= 1
        else
            population.compartment_vars[INFECTED_IMMUNE] += 1
            population.compartment_vars[UNINFECTED_IMMUNE] -= 1
        end
    end

    push!(population.hosts[host_idx].pathogens, pathogen)
    push!(population.hosts[host_idx].pathogen_fractions, 0.0)
    population.hosts[host_idx].pathogen_weights = catCol(
        population.hosts[host_idx].pathogen_weights, zeros(Float64, NUM_PATHOGEN_EVENTS)
    )

    hostWeights!(host_idx, population, model)
end

function addResponseToHost!(response::Response, host_idx::Int64, population::Population, model::Model)
    if length(population.hosts[host_idx].responses) == 0
        if length(population.hosts[host_idx].pathogens) == 0
            population.compartment_vars[UNINFECTED_IMMUNE] += 1
            population.compartment_vars[UNINFECTED_NAIVE] -= 1
        else
            population.compartment_vars[INFECTED_IMMUNE] += 1
            population.compartment_vars[INFECTED_NAIVE] -= 1
        end
    end

    push!(population.hosts[host_idx].responses, response)
    population.hosts[host_idx].response_weights = catCol(
        population.hosts[host_idx].response_weights, zeros(Float64, NUM_RESPONSE_EVENTS)
    )

    hostWeights!(host_idx, population, model)
end

function removePathogenFromHost!(pathogen_idx::Int64, host_idx::Int64, population::Population, model::Model)
    deleteat!(population.hosts[host_idx].pathogens, pathogen_idx)
    deleteat!(population.hosts[host_idx].pathogen_fractions, pathogen_idx)
    population.hosts[host_idx].pathogen_weights = population.hosts[host_idx].pathogen_weights[:, begin:end.!=pathogen_idx]

    if length(population.hosts[host_idx].pathogens) == 0
        if length(population.hosts[host_idx].responses) == 0
            population.compartment_vars[INFECTED_NAIVE] -= 1
            population.compartment_vars[UNINFECTED_NAIVE] += 1
        else
            population.compartment_vars[INFECTED_IMMUNE] -= 1
            population.compartment_vars[UNINFECTED_IMMUNE] += 1
        end
    end

    hostWeights!(host_idx, population, model)
end

function removeResponseFromHost!(response_idx::Int64, host_idx::Int64, population::Population, model::Model)
    deleteat!(population.hosts[host_idx].responses, response_idx)
    population.hosts[host_idx].response_weights = population.hosts[host_idx].response_weights[:, begin:end.!=response_idx]

    if length(population.hosts[host_idx].responses) == 0
        if length(population.hosts[host_idx].pathogens) == 0
            population.compartment_vars[UNINFECTED_IMMUNE] -= 1
            population.compartment_vars[UNINFECTED_NAIVE] += 1
        else
            population.compartment_vars[INFECTED_IMMUNE] -= 1
            population.compartment_vars[INFECTED_NAIVE] += 1
        end
    end

    hostWeights!(host_idx, population, model)
end

function attemptInfection!(pathogen::Pathogen, host_idx::Int64, pop_idx::Int64, model::Model)
    if !(pathogen in model.populations[pop_idx].hosts[host_idx].pathogens) && rand() < nonsamplingValue(
        TRANSMISSION_EFFICIENCY, pathogen,
        model.populations[pop_idx].hosts[host_idx],
        model.populations[pop_idx]
    )

        addPathogenToHost!(
            pathogen, host_idx, model.populations[pop_idx], model
        )
    end
end

function addHostToPopulation!(new_host::Host, population::Population, model::Model)
    if length(new_host.pathogens) == 0
        if length(new_host.responses) == 0
            population.compartment_vars[UNINFECTED_NAIVE] += 1
        else
            population.compartment_vars[UNINFECTED_IMMUNE] += 1
        end
    else
        if length(new_host.responses) == 0
            population.compartment_vars[INFECTED_NAIVE] += 1
        else
            population.compartment_vars[INFECTED_IMMUNE] += 1
        end
    end

    push!(population.hosts, new_host)
    population.host_weights = catCol(population.host_weights, zeros(Float64, NUM_EVENTS))
    population.host_weights_receive = catCol(population.host_weights_receive, zeros(Float64, NUM_CHOICE_MODIFIERS))
    for evt in EVENTS
        push!(population.host_weights_with_coefficient_samplers[evt], 0.0)
    end
    for evt in CHOICE_MODIFIERS
        push!(population.host_weights_receive_with_coefficient_samplers[evt-CHOICE_MODIFIERS[1]+1], 0.0)
    end

    propagateWeightsOnAddHost!(length(population.hosts), population, model)
end

function addHostsToPopulation!(num_hosts::Int64, host_sequence::String, type::HostType, population::Population, model::Model)
    population.compartment_vars[UNINFECTED_NAIVE] += num_hosts

    # update matrices
    population.host_weights = hcat(population.host_weights, zeros(Float64, NUM_EVENTS, num_hosts))
    population.host_weights_receive = hcat(population.host_weights_receive, zeros(Float64, NUM_CHOICE_MODIFIERS, num_hosts))
    for _ in 1:num_hosts
        for evt in EVENTS
            push!(population.host_weights_with_coefficient_samplers[evt], 0.0)
        end
        for evt in CHOICE_MODIFIERS
            push!(population.host_weights_receive_with_coefficient_samplers[evt-CHOICE_MODIFIERS[1]+1], 0.0)
        end
    end

    num_starting_hosts = length(population.hosts)
    host_coefficients = hostSequenceCoefficients(host_sequence, population.id, type)
    for i in 1:num_hosts
        push!(population.hosts, Host(
            length(population.hosts) + 1,
            MVector{2,Union{Host,Nothing}}([nothing, nothing]),
            model.time,
            host_sequence,
            Vector{Pathogen}(undef, 0), Vector{Response}(undef, 0),
            Vector{Float64}(undef, 0),
            Dict{String,MVector{3,Float64}}(),
            Matrix{Float64}(undef, NUM_PATHOGEN_EVENTS, 0),
            Matrix{Float64}(undef, NUM_RESPONSE_EVENTS, 0),
            host_coefficients,
            type,
        ))

        propagateWeightsOnAddHost!(num_starting_hosts + i, population, model)
    end
end

function removeHostFromPopulation!(host_idx::Int64, population::Population, model::Model)
    if length(population.hosts[host_idx].pathogens) == 0
        if length(population.hosts[host_idx].responses) == 0
            population.compartment_vars[UNINFECTED_NAIVE] -= 1
        else
            population.compartment_vars[UNINFECTED_IMMUNE] -= 1
        end
    else
        if length(population.hosts[host_idx].responses) == 0
            population.compartment_vars[INFECTED_NAIVE] -= 1
        else
            population.compartment_vars[INFECTED_IMMUNE] -= 1
        end
    end

    propagateWeightsOnRemoveHost!(host_idx, population, model)

    population.host_weights = @views(population.host_weights[1:end, 1:end.!=host_idx])
    population.host_weights_receive = @views(population.host_weights_receive[1:end, 1:end.!=host_idx])
    for evt in EVENTS
        deleteat!(population.host_weights_with_coefficient_samplers[evt], host_idx)
    end
    for evt in CHOICE_MODIFIERS
        deleteat!(population.host_weights_receive_with_coefficient_samplers[evt-CHOICE_MODIFIERS[1]+1], host_idx)
    end

    deleteat!(population.hosts, host_idx)
end

function setPopulationContactCoefficient!(pop_1::Population, pop_2::Population, coefficient::Float64, model::Model)
    pop_1.population_contact_coefficients[ model.population_dict[pop_2.id] ] = coefficient
    change = (
        pop_1.population_contact_coefficients[ model.population_dict[pop_2.id] ] *
        model.population_weights_receive[ RECEIVE_CONTACT-CHOICE_MODIFIERS[1]+1, model.population_dict[pop_2.id] ] /
        max(
            length(pop_2.hosts) *
            pop_2.parameters.constant_contact_density,
            1.0
        )) - model.population_contact_weights_receive[ model.population_dict[pop_2.id] ,  model.population_dict[pop_1.id] ]
    updatePopulationContactWeightReceiveMatrix!(model.population_dict[pop_1.id], model.population_dict[pop_2.id], change, model)
end

function setPopulationTransitionCoefficient!(pop_1::Population, pop_2::Population, coefficient::Float64, model::Model)
    pop_1.population_transition_coefficients[ model.population_dict[pop_2.id] ] = coefficient
    change = (
        pop_1.population_transition_coefficients[ model.population_dict[pop_2.id] ] *
        model.population_weights_receive[ RECEIVE_TRANSITION-CHOICE_MODIFIERS[1]+1, model.population_dict[pop_2.id] ] /
        max(
            length(pop_2.hosts) *
            pop_2.parameters.constant_transition_density,
            1.0
        )) - model.population_transition_weights_receive[ model.population_dict[pop_2.id], model.population_dict[pop_1.id] ]
    updatePopulationTransitionWeightReceiveMatrix!(model.population_dict[pop_1.id], model.population_dict[pop_2.id], change, model)
end

# Model events

function establishPathogen!(model::Model, rand_n::Float64)
    pathogen_idx, host_idx, pop_idx, rand_n = choosePathogen(PATHOGEN_ESTABLISHMENT, model, rand_n)

    population = model.populations[pop_idx]
    host = population.hosts[host_idx]
    pathogen = population.pathogens[pathogen_idx]

    # We're not multiplying by generations as explained below, but if we were, here's the math:
    # generations = nonsamplingValue(
    #     GENERATIONS_PER_TRANSMISSION, pathogen, host, population
    # ) * nonsamplingValue(
    #     CONTACT, pathogen, host, population
    # ) / nonsamplingValue(
    #     PATHOGEN_ESTABLISHMENT, pathogen, host, population
    # )

    # First decide whether the event was triggered by a mutation or a recombination, then do it
    mutations = nonsamplingValue(MUTATIONS_PER_GENERATION, pathogen, host, population)
    recombinations = nonsamplingValue(RECOMBINATIONS_PER_GENERATION, pathogen, host, population)
    if length(host.pathogens) == 1 || rand_n < (mutations / max(mutations + recombinations, 1.0))
        new_lineage = mutantPathogen!(
            pathogen, host, population, model.time, mutations
            # we don't multiply by generations here assuming the established lineage is
            # the center of a mutant cloud, so we account for only a single generation's worth
            # of mutations at a time
            # We also don't add recombinations for the same reason
        )
        attemptInfection!(new_lineage, host_idx, pop_idx, model)
    else
        # This assumes recombination rates are not symmetrical:
        # once a lineage has decided to recombine, its partner is chosen based on
        # pathogen fraction, not recombination frequency
        p_idx_2, rand_n = choosePathogenFromPopulationFraction(host_idx, pop_idx, model, rand_n)
        # p_idx_2, rand_n = choosePathogen(host_idx_1, pop_idx_1, RECOMBINATIONS_PER_GENERATION, model, rand())
        pathogen_2 = host.pathogens[p_idx_2]
        num_recombinations = (
        ) / 2.0 # assumption: net recombination rate is the average of the two pathogens
        new_lineage = recombinantPathogens!(
            pathogen, host, population, model.time, recombinations,
            nonsamplingValue(
                RECOMBINATIONS_PER_GENERATION, pathogen_2, host, population
            )
            # we don't multiply by generations here assuming the established lineage is
            # the center of a mutant cloud, so we account for only a single generation's worth
            # of recombinations at a time
            # We also don't add mutations for the same reason
        )
        attemptInfection!(new_lineage, host_idx, pop_idx, model)
    end
end

function clearPathogen!(model::Model, rand_n::Float64)
    pathogen_idx, host_idx, pop_idx, rand_n = choosePathogen(CLEARANCE, model, rand_n)

    host = model.populations[pop_idx].hosts[host_idx]
    pathogen = model.populations[pop_idx].hosts[host_idx].pathogens[pathogen_idx]

    if rand_n < host.type_clearance_numbers[pathogen.type.id][1] / (
            host.type_clearance_numbers[pathogen.type.id][2] *
            host.type_clearance_numbers[pathogen.type.id][3]
        )
        num_pathogens = length(host.pathogens)
        for p in 1:num_pathogens
                # check from back to front by indexing with num_pathogens+1-p
                # to avoid indexing errors as length(host.pathogens) changes
            if model.populations[pop_idx].hosts[host_idx].pathogens[num_pathogens+1-p].type == pathogen.type
                removePathogen!(num_pathogens+1-p, host_idx, pop_idx, model, rand())
            end
        end
    else
        removePathogen!(pathogen_idx, host_idx, pop_idx, model, rand())
    end
end

function removePathogen!(pathogen_idx::Int64, host_idx::Int64, pop_idx::Int64, model::Model, rand_n::Float64)
    resp_acq = nonsamplingValue(
        RESPONSE_ACQUISITION_UPON_CLEARANCE,
        model.populations[pop_idx].hosts[host_idx].pathogens[pathogen_idx],
        model.populations[pop_idx].hosts[host_idx],
        model.populations[pop_idx]
    )

    if (resp_acq > 0.0 && rand() < resp_acq)
        acquireResponse!(pathogen_idx, host_idx, pop_idx, model, rand_n)
    end

    removePathogenFromHost!(
        pathogen_idx, host_idx,
        model.populations[pop_idx], model
    )
end

function acquireResponse!(model::Model, rand_n::Float64)
    pathogen_idx, host_idx, pop_idx = choosePathogen(RESPONSE_ACQUISITION, model, rand_n)
    acquireResponse!(pathogen_idx, host_idx, pop_idx, model, rand_n)
end

function acquireResponse!(pathogen_idx::Int64, host_idx::Int64, pop_idx::Int64, model::Model, rand_n::Float64)
    responses = model.populations[pop_idx].parameters.developResponses(
        model.populations[pop_idx].hosts[host_idx].pathogens[pathogen_idx],
        model.populations[pop_idx].hosts[host_idx],
        model.populations[pop_idx].id,
        model.populations[pop_idx].responses,
        model.populations[pop_idx].parameters.response_types,
        model.time
    )

    for response in responses
        if !(response in model.populations[pop_idx].hosts[host_idx].responses)
            addResponseToHost!(
                response, host_idx,
                model.populations[pop_idx], model
            )
        end
    end
end

function hostContact!(model::Model, rand_n::Float64)
    host_idx_1, pop_idx_1, rand_n = chooseHost(CONTACT, model, rand_n)
    pop_idx_2, rand_n = choosePopulationReceiveContact(pop_idx_1, model, rand_n)
    host_idx_2, rand_n = chooseHost(pop_idx_2, RECEIVE_CONTACT, model, rand_n)

    hostContact!(host_idx_1, pop_idx_1, host_idx_2, pop_idx_2, model, rand_n)
end

function hostContact!(
    host_idx_1::Int64, pop_idx_1::Int64, host_idx_2::Int64, pop_idx_2::Int64,
    model::Model, rand_n::Float64; inoculum_coefficient=INOCULUM)
    if host_idx_1 != host_idx_2 || pop_idx_1 != pop_idx_2
        host1 = model.populations[pop_idx_1].hosts[host_idx_1]
        # inocula = MVector{length(model.populations[pop_idx_1].hosts[host_idx_1].pathogens),Int64}([
        inocula = Vector{Int64}([
            pois_rand(
                nonsamplingValue(
                    inoculum_coefficient, host1.pathogens[p_idx], host1,
                    model.populations[pop_idx_1]
                ) * host1.pathogen_fractions[p_idx]
            )
            for p_idx in 1:length(host1.pathogens)
        ])
        if all(inocula .== 0.0) # had originally written !any(inocula .!= 0.0), don't know if that's actually faster
            idx, rand_n = choosePathogenFromPopulationFraction(host_idx_1, pop_idx_1, model, rand_n)
            inocula[idx] += 1
        end

        for p_idx in eachindex(inocula)
            generations = nonsamplingValue(
                GENERATIONS_PER_TRANSMISSION, host1.pathogens[p_idx],
                host1, model.populations[pop_idx_1]
            )
            mutations = nonsamplingValue(
                MUTATIONS_PER_GENERATION, host1.pathogens[p_idx], host1,
                model.populations[pop_idx_1]
            ) * generations
            recombinations = nonsamplingValue(
                RECOMBINATIONS_PER_GENERATION, host1.pathogens[p_idx], host1,
                model.populations[pop_idx_1]
            ) * generations

            if inocula[p_idx] > 0
                mut_prob = min(
                    1.0, 1.0 - exp(
                        -mutations
                    )
                )
                # probability from Poisson PMF with k=0
                rec_prob = 0.0
                if length(host1.pathogens) > 1
                    rec_prob = min(
                        1.0, 1.0 - exp(
                            -recombinations
                        )
                    )
                    # probability from Poisson PMF with k=0
                end
                num_mut = binomial(inocula[p_idx], mut_prob * (1.0 - rec_prob))
                num_rec = binomial(inocula[p_idx], rec_prob * (1.0 - mut_prob))
                num_mut_rec = binomial(inocula[p_idx], mut_prob * rec_prob)
                if inocula[p_idx] - num_mut - num_rec - num_mut_rec > 0
                    attemptInfection!(
                        host1.pathogens[p_idx],
                        host_idx_2, pop_idx_2, model
                    )
                end
                for _ in 1:num_mut
                    attemptInfection!(
                        mutantPathogen!(
                            host1.pathogens[p_idx],
                            host1,
                            model.populations[pop_idx_1],
                            model.time,
                            mutations
                        ), host_idx_2, pop_idx_2, model
                    )
                end
                for _ in 1:num_rec
                    # This assumes recombination rates are not symmetrical:
                    # once a lineage has decided to recombine, its partner is chosen based on
                    # pathogen fraction, not recombination frequency
                    p_idx_2, rand_n = choosePathogenFromPopulationFraction(host_idx_1, pop_idx_1, model, rand_n)
                    # p_idx_2, rand_n = choosePathogen(host_idx_1, pop_idx_1, RECOMBINATIONS_PER_GENERATION, model, rand())

                    if p_idx != p_idx_2 && (
                        host1.pathogens[p_idx].type ==
                        host1.pathogens[p_idx_2].type)

                        attemptInfection!(
                            recombinantPathogens!(
                                host1.pathogens[p_idx],
                                host1.pathogens[p_idx_2],
                                host1,
                                model.populations[pop_idx_1],
                                model.time,
                                recombinations,
                                nonsamplingValue(
                                    RECOMBINATIONS_PER_GENERATION,
                                    host1.pathogens[p_idx_2],
                                    host1, model.populations[pop_idx_1]
                                )
                            ), host_idx_2, pop_idx_2, model
                        )
                    end
                end
                for _ in 1:num_mut_rec
                    # This assumes recombination rates are not symmetrical:
                    # once a lineage has decided to recombine, its partner is chosen based on
                    # pathogen fraction, not recombination frequency
                    p_idx_2, rand_n = p_idx_2, rand_n = choosePathogenFromPopulationFraction(host_idx_1, pop_idx_1, model, rand_n)
                    # p_idx_2, rand_n = choosePathogen(host_idx_1, pop_idx_1, RECOMBINATIONS_PER_GENERATION, model, rand())

                    if p_idx != p_idx_2
                        recombinant = recombinantPathogens!(
                            host1.pathogens[p_idx],
                            host1.pathogens[p_idx_2],
                            host1,
                            model.populations[pop_idx_1],
                            model.time,
                            recombinations,
                            nonsamplingValue(
                                RECOMBINATIONS_PER_GENERATION, host1.pathogens[p_idx_2],
                                host1, model.populations[pop_idx_1]
                            )
                        )

                        attemptInfection!(
                            mutantPathogen!(
                                recombinant, host1,
                                model.populations[pop_idx_1], model.time,
                                mutations
                            ),
                            host_idx_2, pop_idx_2, model
                        )
                    else
                        attemptInfection!(
                            mutantPathogen!(
                                host1.pathogens[p_idx],
                                host1,
                                model.populations[pop_idx_1],
                                model.time,
                                mutations
                            ),
                            host_idx_2, pop_idx_2, model
                        )
                    end
                end
            end
        end
    end
end

function loseResponse!(model::Model, rand_n::Float64)
    response_idx, host_idx, pop_idx = chooseResponse(RESPONSE_LOSS, model, rand_n)

    removeResponseFromHost!(
        response_idx, host_idx, model.populations[pop_idx], model
    )
end

function birth!(model::Model, rand_n::Float64)
    host_idx, pop_idx, rand_n = chooseHost(BIRTH, model, rand_n)

    parents = MVector{2,Union{Nothing,Host}}(model.populations[pop_idx].hosts[host_idx], nothing)
    if model.populations[pop_idx].parameters.host_sexual_reproduction
        host_idx_2, pop_idx_2, rand_n = chooseHost(BIRTH, model, rand_n)
        if model.populations[pop_idx].parameters.hostSexualCompatibility(parents[1].sequence, model.populations[pop_idx_2].hosts[host_idx_2].sequence)
            parents = MVector{2,Union{Nothing,Host}}(
                model.populations[pop_idx].hosts[host_idx],
                model.populations[pop_idx_2].hosts[host_idx_2]
            )
        else
            parents = MVector{2,Union{Nothing,Host}}(nothing, nothing)
        end
    end

    if !isnothing(parents[1])
        child_sequence = ""

        parent_gametes = MVector{2,Union{String,Nothing}}(
            generateGamete(
                parents[1].sequence,
                parents[1].type.num_loci,
                parents[1].type.possible_alleles,
                nonsamplingValue(
                    HOST_MUTATIONS_UPON_BIRTH,
                    parents[1],
                    model.populations[pop_idx]
                )
            ), nothing
        )

        if occursin(HOMOLOGOUS_CHROMOSOME_SEPARATOR, parent_gametes[1]) # if diploid, #model.populations[pop_idx].diploid_host
            parent_gametes[1] = recombinantSequences(
                parent_gametes[1],
                generateGamete(
                    parents[1].sequence,
                    parents[1].type.num_loci,
                    parents[1].type.possible_alleles,
                    nonsamplingValue(
                        HOST_MUTATIONS_UPON_BIRTH,
                        parents[1],
                        model.populations[pop_idx]
                    )
                ),
                parents[1].type.num_loci,
                parents[1].mean_recombination_crossovers, #TODO: this parameter doesn't exist anymore?
                parents[1].mean_recombination_crossovers, #TODO: this parameter doesn't exist anymore?
            )
            if !isnothing(parents[2])
                parent_gametes[2] = recombinantSequences(
                    parent_gametes[2],
                    generateGamete(
                        parents[2].sequence,
                        parents[2].type.num_loci,
                        parents[2].type.possible_alleles,
                        nonsamplingValue(
                            HOST_MUTATIONS_UPON_BIRTH,
                            parents[2],
                            model.populations[pop_idx]
                        )
                    ),
                    parents[2].type.num_loci,
                    parents[2].mean_recombination_crossovers, #TODO: this parameter doesn't exist anymore?
                    parents[2].mean_recombination_crossovers #TODO: this parameter doesn't exist anymore?
                )
            end
            child_sequence = generateZygote(parent_gametes[1], parent_gametes[2])
        else
            rec_1 = nonsamplingValue(
                HOST_RECOMBINATIONS_UPON_BIRTH,
                parents[1],
                model.populations[pop_idx]
            )
            rec_2 = nonsamplingValue(
                HOST_RECOMBINATIONS_UPON_BIRTH,
                parents[2],
                model.populations[pop_idx]
            )
            if !isnothing(parents[2]) && rec_1 > 0.0 && rec_2 > 0.0
                child_sequence = recombinantSequences(
                    parent_gametes[1],
                    generateGamete(
                        parents[2].sequence,
                        parents[2].type.num_loci,
                        parents[2].type.possible_alleles,
                        nonsamplingValue(
                            HOST_MUTATIONS_UPON_BIRTH,
                            parents[2],
                            model.populations[pop_idx]
                        )
                    ),
                    parents[2].type.num_loci,
                    rec_1, rec_2
                )
            else
                child_sequence = parent_gametes[1]
            end
        end

        jOpqua.newHost!(
            child_sequence, parents[1].type, model.populations[pop_idx], model,
            parents=parents, birth_time=model.time
        )

        for parent in parents
            if !isnothing(parent)
                for response in parent.responses
                    resp_inh = nonsamplingValue(
                        RESPONSE_INHERITANCE, response, parent, model.populations[pop_idx]
                    )
                    if resp_inh > 0.0 && rand() < resp_inh
                        addResponseToHost!(
                            response, length(model.populations[pop_idx].hosts),
                            model.populations[pop_idx], model
                        )
                    end
                end

                vert_trans = nonsamplingValue(
                    VERTICAL_TRANSMISSION, parent,
                    model.populations[pop_idx]
                )

                if (vert_trans > 0.0 && rand() < vert_trans)
                    hostContact!(
                        parent, pop_idx, length(model.populations[pop_idx].hosts),
                        pop_idx, model, rand_n, inoculum_coefficient=VERTICAL_TRANSMISSION
                    )
                end
            end
        end
    end
end

function death!(model::Model, rand_n::Float64)
    host_idx, pop_idx, rand_n = chooseHost(DEATH, model, rand_n)

    model.populations[pop_idx].compartment_vars[DEAD] += 1
    # done here because removeHostFromPopulation doesn't advance deaths,
    # but does remove from other vars

    removeHostFromPopulation!(host_idx, model.populations[pop_idx], model)
end

function transition!(model::Model, rand_n::Float64)
    host_idx, pop_idx_1, rand_n = chooseHost(TRANSITION, model, rand_n)
    pop_idx_2, rand_n = choosePopulationReceiveTransition(pop_idx_1, model, rand_n)

    addHostToPopulation!(model.populations[pop_idx_1].hosts[host_idx], model.populations[pop_idx_2], model)
    removeHostFromPopulation!(host_idx, model.populations[pop_idx_1], model)
end

const EVENT_FUNCTIONS = SA[
    establishPathogen!, clearPathogen!, acquireResponse!,
    hostContact!, loseResponse!,
    birth!, death!, transition!
]
