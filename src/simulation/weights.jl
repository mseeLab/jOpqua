using StaticArrays
using Flexle

## Coefficient and weight calculation (from bottom to top):

# Bottom-level coefficients (Response and Pathogen)

function responseStaticSpecificCoefficients(
    imprinted_pathogen_sequence::String, matured_pathogen_sequence::String,
    host_sequence::String, population_id::String, type::ResponseType)
    return [
        type.static_specific_coefficient_functions[evt_id](
            imprinted_pathogen_sequence, matured_pathogen_sequence, host_sequence, population_id
        )
        for evt_id in 1:NUM_COEFFICIENTS
    ]
end

function responseStaticHostwideCoefficients(
    imprinted_pathogen_sequence::String, matured_pathogen_sequence::String,
    host_sequence::String, population_id::String, type::ResponseType)
    return [
        type.static_hostwide_coefficient_functions[evt_id](
            imprinted_pathogen_sequence, matured_pathogen_sequence, host_sequence, population_id
        )
        for evt_id in 1:NUM_COEFFICIENTS
    ]
end

function responseStaticSpecificCoefficient(response::Response, host::Host, population::Population, coefficient::Int64)
    return response.type.static_specific_coefficient_functions[coefficient](
        isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
        isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
        host.sequence,
        population.id
    )
end

function responseInteractionSpecificCoefficient(pathogen::Pathogen, response::Response, host::Host, population::Population, coefficient::Int64)
    return response.type.interaction_specific_coefficient_functions[coefficient](
        pathogen.sequence,
        isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
        isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
        host.sequence,
        population.id
    )
end

function responseStaticHostwideCoefficient(response::Response, host::Host, population::Population, coefficient::Int64)
    return response.type.static_hostwide_coefficient_functions[coefficient](
        isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
        isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
        host.sequence,
        population.id
    )
end

function responseInteractionHostwideCoefficient(pathogen::Pathogen, response::Response, host::Host, population::Population, coefficient::Int64)
    return response.type.interaction_hostwide_coefficient_functions[coefficient](
        pathogen.sequence,
        isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
        isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
        host.sequence,
        population.id
    )
end

function reactivityCoefficient(pathogen::Pathogen, response::Response, host::Host, population::Population)
    return response.type.reactivityCoefficient(
        pathogen.sequence,
        isnothing(response.imprinted_pathogen) ? "" : response.imprinted_pathogen.sequence,
        isnothing(response.matured_pathogen) ? "" : response.matured_pathogen.sequence,
        host.sequence,
        population.id
    )
end

function nonsamplingValue(coef::Int64, host::Host, population::Population)
    return host.coefficients[coef] * hostwideNetCoefficient(host, population, coef) * population.parameters.base_coefficients[coef]
end

function nonsamplingValue(coef::Int64, pathogen::Pathogen, host::Host, population::Population)
    return pathogen.specific_coefficients[coef] * population.parameters.weightedInteractionPathogen(
               pathogen, host, coef
           ) * host.coefficients[coef] * hostwideNetCoefficient(host, population, coef) * population.parameters.base_coefficients[coef]
end

function nonsamplingValue(coef::Int64, response::Response, host::Host, population::Population)
    return response.specific_coefficients[coef] * population.parameters.weightedInteractionResponse(
               response, host, coef
           ) * host.coefficients[coef] * hostwideNetCoefficient(host, population, coef) * population.parameters.base_coefficients[coef]
end

function pathogenSequenceSpecificCoefficients(sequence::String, population_id::String, type::PathogenType)
    return [
        type.specific_coefficient_functions[evt_id](sequence, population_id)
        for evt_id in 1:NUM_COEFFICIENTS
    ]
end

function pathogenSequenceHostwideCoefficients(sequence::String, population_id::String, type::PathogenType)
    return [
        type.hostwide_coefficient_functions[evt_id](sequence, population_id)
        for evt_id in 1:NUM_COEFFICIENTS
    ]
end

function hostSequenceCoefficients(sequence::String, population_id::String, type::HostType)
    return MVector{NUM_COEFFICIENTS,Float64}([
        type.coefficient_functions[evt_id](sequence, population_id)
        for evt_id in 1:NUM_COEFFICIENTS
    ])
end

# Intrahost-level weights

function hostwideNetCoefficient(host::Host, population::Population, coef::Int64)
    net_hostwide_coefficient = 1.0
    if length(host.pathogens) > 0
        net_hostwide_coefficient = 0.0
        for p in 1:length(host.pathogens)
            net_hostwide_coefficient += host.pathogens[p].hostwide_coefficients[coef] * host.pathogen_fractions[p]
        end
    end
    for re in 1:length(host.responses)
        net_hostwide_coefficient = net_hostwide_coefficient * host.responses[re].hostwide_coefficients[coef]
    end

    return net_hostwide_coefficient * population.parameters.weightedInteractionHostwide(host, coef) * host.coefficients[coef]
end

function pathogenWeights!(p::Int64, host::Host, population::Population, evt::Int64)
    host.pathogen_weights[evt, p] = host.pathogens[p].specific_coefficients[evt]
    if length(host.responses) > 0
        host.pathogen_weights[evt, p] =
            host.pathogen_weights[evt, p] * population.parameters.weightedInteractionPathogen(
                host.pathogens[p], host, evt
            )
    end
    if evt in PATHOGEN_NONSAMPLING_COEFFICIENTS
        # PATHOGEN_NONSAMPLING_COEFFICIENTS don't require host-level aggregation
    elseif evt == CLEARANCE
        # Clearance likelihoods are not proportional to intrahost fraction,
        # instead, we assume the rate of loss compounds
        # (could be inversely proportional to fraction?)
        if !haskey(host.type_clearance_numbers, host.pathogens[p].type.id)
            host.type_clearance_numbers[host.pathogens[p].type.id] = MVector{3,Float64}(zeros(Float64, 3))
        end
        host.type_clearance_numbers[host.pathogens[p].type.id][1] += host.pathogen_weights[evt, p]
        host.type_clearance_numbers[host.pathogens[p].type.id][2] += 1.0
        host.pathogen_weights[evt, p] += pathogenClearanceThroughDriftWeight(host.pathogens[p], host, population)
        host.type_clearance_numbers[host.pathogens[p].type.id][3] += host.pathogen_weights[evt, p]
    elseif evt == PATHOGEN_ESTABLISHMENT && length(host.pathogens) == 1
        # if nobody to recombine with, no recombination lineages can establish, so we
        # subtract the fraction of the weight that came from recombinations
        host.pathogen_weights[evt, p] =
            host.pathogen_weights[evt, p] *
            host.pathogen_fractions[p] -
            nonsamplingValue(RECOMBINATIONS_PER_GENERATION, host.pathogens[p], host, population) / max(
                nonsamplingValue(MUTATIONS_PER_GENERATION, host.pathogens[p], host, population) +
                nonsamplingValue(RECOMBINATIONS_PER_GENERATION, host.pathogens[p], host, population), 1.0
            )
        #TODO: regarding recombination and mutation establishment, each of these
        # depend on other parameters within Pathogen (mean_mutations_per_replication
        # and mean_recombination_crossovers); since in the current version we don't
        # consider intrahost population sizes or replication rates, then it's
        # pointless to incorporate these parameters in weight calculation here, but
        # eventually in a future version with better intrahost dynamics and popgen
        # this will be the case
    else
        # All other pathogen event likelihoods are proportional to intrahost fraction
        host.pathogen_weights[evt, p] =
            host.pathogen_weights[evt, p] *
            host.pathogen_fractions[p]
    end
end

function hostWeightsPathogen!(host_idx::Int64, population::Population, evt::Int64)
    population.host_weights[evt, host_idx] = 0.0
    for p in 1:length(population.hosts[host_idx].pathogens)
        pathogenWeights!(p, population.hosts[host_idx], population, evt)
        population.host_weights[evt, host_idx] += population.hosts[host_idx].pathogen_weights[evt, p]
        # we sum here because we have already weighted by fraction
        # (for everything except clearance, see above)
    end
    # if evt == CLEARANCE
    #     for p in 1:length(population.hosts[host_idx].pathogens)
    #         pathogenWeights!(p, population.hosts[host_idx], population, evt)
    #         if population.hosts[host_idx].pathogen_weights[evt, p] > population.host_weights[evt, host_idx]
    #             population.host_weights[evt, host_idx] = population.hosts[host_idx].pathogen_weights[evt, p]
    #         end
    #         # For clearance, we take the maximum rate
    #     end
    # else
    #     for p in 1:length(population.hosts[host_idx].pathogens)
    #         pathogenWeights!(p, population.hosts[host_idx], population, evt)
    #         population.host_weights[evt, host_idx] += population.hosts[host_idx].pathogen_weights[evt, p]
    #         # we sum here because we have already weighted by fraction
    #         # (for everything except clearance, see above)
    #     end
    # end

    population.host_weights[evt, host_idx] = (
        population.host_weights[evt, host_idx] *
        hostwideNetCoefficient(population.hosts[host_idx], population, evt)
    )

    setindex!(
        population.host_weights_with_coefficient_samplers[evt],
        population.host_weights[evt, host_idx] *
        population.parameters.base_coefficients[evt],
        host_idx
    )
end

function responseWeights!(re::Int64, host::Host, population::Population, evt::Int64)
    host.response_weights[evt-RESPONSE_EVENTS[1]+1, re] = host.responses[re].specific_coefficients[evt] * population.parameters.weightedInteractionResponse(
        host.responses[re], host, evt
    )
    # No Response fraction weighting here, just presence of Response is enough
end

function hostWeightsResponse!(host_idx::Int64, population::Population, evt::Int64)
    population.host_weights[evt, host_idx] = 0.0
    for re in 1:length(population.hosts[host_idx].responses)
        responseWeights!(re, population.hosts[host_idx], population, evt)
        population.host_weights[evt, host_idx] += population.hosts[host_idx].response_weights[evt-RESPONSE_EVENTS[1]+1, re]
        # We sum here because having more responses does imply more events, e.g. losses of response
    end

    population.host_weights[evt, host_idx] = (
        population.host_weights[evt, host_idx] *
        hostwideNetCoefficient(population.hosts[host_idx], population, evt)
    )

    setindex!(
        population.host_weights_with_coefficient_samplers[evt],
        population.host_weights[evt, host_idx] *
        population.parameters.base_coefficients[evt],
        host_idx
    )
end

# Intra-Population level

function hostWeightsHost!(h::Int64, population::Population, evt::Int64)
    population.host_weights[evt, h] = 1.0 * hostwideNetCoefficient(population.hosts[h], population, evt)
    setindex!(
        population.host_weights_with_coefficient_samplers[evt],
        population.host_weights[evt, h] *
        population.parameters.base_coefficients[evt],
        h
    )
end

function hostWeightsReceive!(h::Int64, population::Population, evt::Int64)
    population.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] = 1.0 * hostwideNetCoefficient(population.hosts[h], population, evt)
    setindex!(
        population.host_weights_receive_with_coefficient_samplers[evt-CHOICE_MODIFIERS[1]+1],
        population.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, h] *
        population.parameters.base_coefficients[evt],
        h
    )
end

function hostWeights!(host_idx::Int64, population::Population, model::Model)
    population.hosts[host_idx].pathogen_fractions = population.parameters.pathogenFractions(
        population.hosts[host_idx], population.parameters.weightedInteractionPathogen
    )
    prev = 0.0
    for weight in PATHOGEN_EVENTS
        prev = population.host_weights[weight, host_idx]
        hostWeightsPathogen!(host_idx, population, weight)
        if prev != population.host_weights[weight, host_idx] && population.parameters.base_coefficients[weight] != 0.0
            change = population.host_weights[weight, host_idx] - prev
            if weight == CONTACT
                population.contact_sum += change
                change = change * model.population_contact_weights_receive_sums[model.population_dict[population.id]]
            end
            propagateWeightChanges!(
                population.parameters.base_coefficients[weight] *
                change,
                population, weight, model
            )
        end
    end
    for weight in RESPONSE_EVENTS
        prev = population.host_weights[weight, host_idx]
        hostWeightsResponse!(host_idx, population, weight)
        if prev != population.host_weights[weight, host_idx] && population.parameters.base_coefficients[weight] != 0.0
            propagateWeightChanges!(
                population.parameters.base_coefficients[weight] *
                (population.host_weights[weight, host_idx] - prev),
                population, weight, model
            )
        end
    end
    for weight in HOST_EVENTS
        prev = population.host_weights[weight, host_idx]
        hostWeightsHost!(host_idx, population, weight)
        if prev != population.host_weights[weight, host_idx] && population.parameters.base_coefficients[weight] != 0.0
            change = population.host_weights[weight, host_idx] - prev
            if weight == TRANSITION
                population.transition_sum += change
                change = change * model.population_transition_weights_receive_sums[model.population_dict[population.id]]
            end
            propagateWeightChanges!(
                population.parameters.base_coefficients[weight] *
                change,
                population, weight, model
            )
        end
    end
    for weight in CHOICE_MODIFIERS[begin:end-1]
        prev = population.host_weights_receive[weight-CHOICE_MODIFIERS[1]+1, host_idx]
        hostWeightsReceive!(host_idx, population, weight)
        if (
            prev != population.host_weights_receive[weight-CHOICE_MODIFIERS[1]+1, host_idx] &&
            population.parameters.base_coefficients[weight] != 0.0 &&
            weight < INTRAHOST_FITNESS
        )
            propagateWeightReceiveChanges!(
                population.parameters.base_coefficients[weight] *
                (population.host_weights_receive[weight-CHOICE_MODIFIERS[1]+1, host_idx] - prev),
                population, weight, model
            )
        end
    end
end

# Intra-Event/Weight level:

function calculateWeight!(population::Population, evt::Int64, model::Model)
    if evt == CONTACT
        population.contact_sum = sum(@views(population.host_weights[evt, :]))
        # The latest version of Flexle performs the above error correction internally
        model.population_weights[evt, model.population_dict[population.id]] = (
            model.population_contact_weights_receive_sums[model.population_dict[population.id]] *
            population.parameters.base_coefficients[evt] * population.contact_sum
        )
    elseif evt == TRANSITION
        population.transition_sum = sum(@views(population.host_weights[evt, :]))
        model.population_weights[evt, model.population_dict[population.id]] = (
            model.population_transition_weights_receive_sums[model.population_dict[population.id]] *
            population.parameters.base_coefficients[evt] * population.transition_sum
        )
        # The latest version of Flexle performs the above error correction internally
    else
        model.population_weights[evt, model.population_dict[population.id]] = population.parameters.base_coefficients[evt] * sum(@views(population.host_weights[evt, :]))
        # The latest version of Flexle performs the above error correction internally
    end
    population.recalculation_counters[evt] = 0
end

function calculateWeightReceive!(population::Population, weight::Int64, model::Model)
    if length(population.hosts) == 0 && weight == RECEIVE_TRANSITION
        model.population_weights_receive[weight-CHOICE_MODIFIERS[1]+1, model.population_dict[population.id]] = population.parameters.base_coefficients[weight]
        # receiving transitions doesn't stop happening with an empty population
    else
        model.population_weights_receive[weight-CHOICE_MODIFIERS[1]+1, model.population_dict[population.id]] = population.parameters.base_coefficients[weight] * sum(@views(population.host_weights_receive[weight-CHOICE_MODIFIERS[1]+1, :]))
    end
    population.recalculation_counters[weight] = 0
end

function calculatePopulationContactWeightReceiveMatrix!(pop_idx_i::Int64, model::Model)
    model.population_contact_weights_receive_sums[pop_idx_i] = 0
    for pop_idx_j in values(model.population_dict)
        model.population_contact_weights_receive[pop_idx_j, pop_idx_i] = (
            model.populations[pop_idx_i].population_contact_coefficients[pop_idx_j] *
            model.population_weights_receive[RECEIVE_CONTACT-CHOICE_MODIFIERS[1]+1, pop_idx_j] /
            max(
                length(model.populations[pop_idx_j].hosts) *
                model.populations[pop_idx_j].parameters.constant_contact_density,
                1.0
            )
        )
        model.population_contact_weights_receive_sums[pop_idx_i] += model.population_contact_weights_receive[pop_idx_j, pop_idx_i]
    end
end

function calculatePopulationTransitionWeightReceiveMatrix!(pop_idx_i::Int64, model::Model)
    model.population_transition_weights_receive_sums[pop_idx_i] = 0
    for pop_idx_j in values(model.population_dict)
        model.population_transition_weights_receive[pop_idx_j, pop_idx_i] = (
            model.populations[pop_idx_i].population_transition_coefficients[pop_idx_j] *
            model.population_weights_receive[RECEIVE_TRANSITION-CHOICE_MODIFIERS[1]+1, pop_idx_j] /
            max(
                length(model.populations[pop_idx_j].hosts) *
                model.populations[pop_idx_j].parameters.constant_transition_density,
                1.0
            )
        )
        model.population_transition_weights_receive_sums[pop_idx_i] += model.population_transition_weights_receive[pop_idx_j, pop_idx_i]
    end
end

function updatePopulationContactWeightReceiveMatrix!(pop_idx_1::Int64, pop_idx_2::Int64, change::Float64, model::Model)
    model.population_contact_weights_receive[pop_idx_2, pop_idx_1] += change
    if becomesZero(model.population_contact_weights_receive_sums[pop_idx_1], change)
        prev = model.population_contact_weights_receive_sums[pop_idx_1]
        calculatePopulationContactWeightReceiveMatrix!(pop_idx_1, model)
        change = model.population_contact_weights_receive_sums[pop_idx_1] - prev
    end

    propagateWeightChanges!(
        model.populations[pop_idx_1].parameters.base_coefficients[CONTACT] *
        model.populations[pop_idx_1].contact_sum * change,
        model.populations[pop_idx_1], CONTACT, model
    )
end

function updatePopulationTransitionWeightReceiveMatrix!(pop_idx_1::Int64, pop_idx_2::Int64, change::Float64, model::Model)
    model.population_transition_weights_receive[pop_idx_2, pop_idx_1] += change
    if becomesZero(model.population_transition_weights_receive_sums[pop_idx_1], change)
        prev = model.population_transition_weights_receive_sums[pop_idx_1]
        calculatePopulationTransitionWeightReceiveMatrix!(pop_idx_1, model)
        change = model.population_transition_weights_receive_sums[pop_idx_1] - prev
    end

    propagateWeightChanges!(
        model.populations[pop_idx_1].parameters.base_coefficients[TRANSITION] *
        model.populations[pop_idx_1].transition_sum * change,
        model.populations[pop_idx_1], TRANSITION, model
    )
end

# Intra-model level:

function calculateRate!(evt::Int64, model::Model)
    model.event_rates[evt] = sum([
        model.population_weights[evt, pop_idx] for pop_idx in values(model.population_dict)
    ])
    model.recalculation_counters[evt] = 0
end

function calculateWeightReceive!(weight::Int64, model::Model)
    model.population_weights_receive_sums[weight-CHOICE_MODIFIERS[1]+1] = sum([
        model.population_weights_receive[weight-CHOICE_MODIFIERS[1]+1, pop_idx] for pop_idx in values(model.population_dict)
    ])
    model.recalculation_counters[weight] = 0
end

# Weight change propagation:

function propagateWeightChanges!(change::Float64, evt::Int64, model::Model)
    if becomesZero(model.event_rates[evt], change) || model.recalculation_counters[evt] > RECALCULATION_EVENT_THRESHOLD
        calculateRate!(evt, model)
    else
        model.event_rates[evt] += change
        model.recalculation_counters[evt] += 1
    end
    model.event_rates_sum = sum(model.event_rates) # this single line is all the "super-model level" calculation code haha
end

function propagateWeightChanges!(change::Float64, population::Population, evt::Int64, model::Model)
    if change != 0.0
        population.recalculation_counters[evt] += 1
        if population.recalculation_counters[evt] < RECALCULATION_EVENT_THRESHOLD
            if becomesZero(model.population_weights[evt, model.population_dict[population.id]], change)
                prev = model.population_weights[evt, model.population_dict[population.id]]
                calculateWeight!(population, evt, model)
                change = model.population_weights[evt, model.population_dict[population.id]] - prev
            else
                model.population_weights[evt, model.population_dict[population.id]] += change
            end
            propagateWeightChanges!(change, evt, model)
        else
            calculateWeight!(population, evt, model)
            calculateRate!(evt, model)
        end
    end
end

function propagateWeightChanges!(change::Float64, host_idx::Int64, population::Population, evt::Int64, model::Model)
    if change != 0
        population.host_weights[evt, host_idx] += change

        if length(population.hosts) > 0
            if evt == CONTACT
                population.contact_sum += change
                change = change * model.population_contact_weights_receive_sums[model.population_dict[population.id]]
            elseif evt == TRANSITION
                population.transition_sum += change
                change = change * model.population_transition_weights_receive_sums[model.population_dict[population.id]]
            end
        else
            if evt == CONTACT
                population.contact_sum = 0.0
            elseif evt == TRANSITION
                population.transition_sum = 0.0
            end
            change = -model.population_weights[evt, model.population_dict[population.id]]
        end

        if population.parameters.base_coefficients[evt] != 0
            setindex!(
                population.host_weights_with_coefficient_samplers[evt],
                population.host_weights[evt, host_idx] *
                population.parameters.base_coefficients[evt],
                host_idx
            )

            propagateWeightChanges!(population.parameters.base_coefficients[evt] * change, population, evt, model)
        end
    end
end

function propagateWeightReceiveChanges!(change::Float64, population::Population, evt::Int64, model::Model)
    if change != 0.0
        population.recalculation_counters[evt] += 1
        model.recalculation_counters[evt] += 1
        if population.recalculation_counters[evt] < RECALCULATION_EVENT_THRESHOLD
            if becomesZero(model.population_weights_receive[evt-CHOICE_MODIFIERS[1]+1, model.population_dict[population.id]], change)
                prev = model.population_weights_receive[evt-CHOICE_MODIFIERS[1]+1, model.population_dict[population.id]]
                calculateWeightReceive!(population, evt, model)
                change = model.population_weights_receive[evt-CHOICE_MODIFIERS[1]+1, model.population_dict[population.id]] - prev
            else
                model.population_weights_receive[evt-CHOICE_MODIFIERS[1]+1, model.population_dict[population.id]] += change
            end
            if ( becomesZero(model.population_weights_receive_sums[evt-CHOICE_MODIFIERS[1]+1], change)
                || model.recalculation_counters[evt] > RECALCULATION_EVENT_THRESHOLD )
                calculateWeightReceive!(evt, model)
            else
                model.population_weights_receive_sums[evt-CHOICE_MODIFIERS[1]+1] += change
            end

            if evt == RECEIVE_CONTACT
                for p in 1:length(model.populations)
                    change_p = (
                        change * model.populations[p].population_contact_coefficients[model.population_dict[population.id]] /
                        max(length(population.hosts) * population.parameters.constant_contact_density, 1.0)
                    )
                    # Contact rates assume scaling area if constant_contact_density is true:
                    # large populations are equally
                    # dense as small ones, so contact is constant (divide by total hosts).
                    # If you don't want this to happen, modify each population's
                    # receive contact coefficient accordingly.
                    model.population_contact_weights_receive[model.population_dict[population.id], p] += change_p
                    model.population_contact_weights_receive_sums[p] += change_p
                    propagateWeightChanges!(
                        change_p * model.populations[p].parameters.base_coefficients[CONTACT] * model.populations[p].contact_sum,
                        model.populations[p], CONTACT, model
                    )
                end
            elseif evt == RECEIVE_TRANSITION
                for p in 1:length(model.populations)
                    change_p = (
                        change * model.populations[p].population_transition_coefficients[model.population_dict[population.id]] /
                        max(length(population.hosts) * population.parameters.constant_transition_density, 1.0)
                    )
                    # Transition receive weights are assumed to be independent of population
                    # (constant_transition_density=true), but can be modified if desired.
                    model.population_transition_weights_receive[model.population_dict[population.id], p] += change_p
                    model.population_transition_weights_receive_sums[p] += change_p
                    propagateWeightChanges!(
                        change_p * model.populations[p].parameters.base_coefficients[TRANSITION] * model.populations[p].transition_sum,
                        model.populations[p], TRANSITION, model
                    )
                end
            end
        else
            calculateWeightReceive!(population, evt, model)
            calculateWeightReceive!(evt, model)

            if evt == RECEIVE_CONTACT
                for p in 1:length(model.populations)
                    calculatePopulationContactWeightReceiveMatrix!(p, model)
                    calculateWeight!(population, CONTACT, model)
                end
            elseif evt == RECEIVE_TRANSITION
                for p in 1:length(model.populations)
                    calculatePopulationTransitionWeightReceiveMatrix!(p, model)
                    calculateWeight!(population, TRANSITION, model)
                end
            end
        end
    end
end

function propagateWeightReceiveChanges!(change::Float64, host_idx::Int64, population::Population, evt::Int64, model::Model)
    if change != 0.0
        population.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, host_idx] += change
        setindex!(
            population.host_weights_receive_with_coefficient_samplers[evt-CHOICE_MODIFIERS[1]+1],
            population.host_weights_receive[evt-CHOICE_MODIFIERS[1]+1, host_idx] *
            population.parameters.base_coefficients[evt],
            host_idx
        )

        if evt < INTRAHOST_FITNESS && population.parameters.base_coefficients[evt] != 0.0
            propagateWeightReceiveChanges!(population.parameters.base_coefficients[evt] * change, population, evt, model)
        end
    end
end

function propagateWeightsOnAddHost!(host_idx::Int64, population::Population, model::Model)
    for p in 1:length(model.populations)
        # First change denominator of contact and transition rate, if necessary
        if population.parameters.constant_contact_density
            change = model.population_contact_weights_receive[model.population_dict[population.id], p] / max(host_idx, 1.0)
            model.population_contact_weights_receive[model.population_dict[population.id], p] -= change
            model.population_contact_weights_receive_sums[p] -= change
            propagateWeightChanges!(
                -change * model.populations[p].parameters.base_coefficients[CONTACT] * model.populations[p].contact_sum,
                model.populations[p], CONTACT, model
            )
        end
        if population.parameters.constant_transition_density || length(population.hosts) == 1
                # if population was empty, receiving transitions are still possible,
                # so make everything zero first before adding new weight later below
            if length(population.hosts) == 1
                model.population_weights_receive[RECEIVE_TRANSITION-CHOICE_MODIFIERS[1]+1, model.population_dict[population.id]] = 0.0
            end
            change = model.population_transition_weights_receive[model.population_dict[population.id], p] / max(host_idx, 1.0)
            model.population_transition_weights_receive[model.population_dict[population.id], p] -= change
            model.population_transition_weights_receive_sums[p] -= change

            propagateWeightChanges!(
                -change * model.populations[p].parameters.base_coefficients[TRANSITION] * model.populations[p].transition_sum,
                model.populations[p], TRANSITION, model
            )
        end
    end

    # Then update the numerator of the rates as if it were a change in the weight of a preexisting host
    hostWeights!(length(population.hosts), population, model)
end

function propagateWeightsOnRemoveHost!(host_idx::Int, population::Population, model::Model)
    # First update the numerator of the rates as if it were a change in the weight of a normal, existing host
    for coef in EVENTS
        if population.host_weights[coef, host_idx] != 0.0
            change = -population.host_weights[coef, host_idx]
            if coef == CONTACT
                population.contact_sum += change
                change = change * model.population_contact_weights_receive_sums[model.population_dict[population.id]]
            elseif coef == TRANSITION
                population.transition_sum += change
                change = change * model.population_transition_weights_receive_sums[model.population_dict[population.id]]
            end
            propagateWeightChanges!(
                change * population.parameters.base_coefficients[coef], population, coef, model
            )
        end
    end
    for coef in CHOICE_MODIFIERS[begin:end-1]
        if population.host_weights_receive[coef-CHOICE_MODIFIERS[1]+1, host_idx] != 0.0
            propagateWeightReceiveChanges!(
                -population.host_weights_receive[coef-CHOICE_MODIFIERS[1]+1, host_idx] * population.parameters.base_coefficients[coef],
                population, coef, model
            )
        end
    end

    # If population is going to be empty, the population receive transition rate must be set to the base receive transition rate
    if length(population.hosts) == 1
        propagateWeightReceiveChanges!(
            population.parameters.base_coefficients[RECEIVE_TRANSITION], population, RECEIVE_TRANSITION, model
        )
    end

    # Then change denominator of contact and transition rate, if necessary
    for p in 1:length(model.populations)
        if population.parameters.constant_contact_density
            change = model.population_contact_weights_receive[model.population_dict[population.id], p] / max(length(population.hosts), 1)
            model.population_contact_weights_receive[model.population_dict[population.id], p] += change
            model.population_contact_weights_receive_sums[p] += change
            propagateWeightChanges!(
                change * model.populations[p].parameters.base_coefficients[CONTACT] * model.populations[p].contact_sum,
                model.populations[p], CONTACT, model
            )
        end
        if population.parameters.constant_transition_density && length(population.hosts) > 1
                # if population is going to be empty, no adjustment of denominator is needed
            change = model.population_transition_weights_receive[model.population_dict[population.id], p] / max(length(population.hosts), 1)
            model.population_transition_weights_receive[model.population_dict[population.id], p] += change
            model.population_transition_weights_receive_sums[p] += change
            propagateWeightChanges!(
                change * model.populations[p].parameters.base_coefficients[TRANSITION] * model.populations[p].transition_sum,
                model.populations[p], TRANSITION, model
            )
        end
    end
end
