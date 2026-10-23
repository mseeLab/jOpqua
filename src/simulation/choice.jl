using Flexle

# Generic sampling methods:
function randChoose(rand_n::Float64, rates::AbstractVector{Float64}, rates_sum::Float64; regenerate_rand::Bool=false)
    #TODO: rates is not type-defined to accommodate @views, maybe this is bad
    if rates_sum > 0.0
        norm_rand_n = rand_n * rates_sum
        cum_sum = 0.0
        idx = 0
        chosen_rate = 0.0
        for r in rates
            idx += 1
            cum_sum += r
            chosen_rate = r
            if cum_sum > norm_rand_n
                break
            end
        end

        if regenerate_rand
            return idx, (norm_rand_n - cum_sum + chosen_rate) / chosen_rate
        else
            return idx, 0.0
        end
    else
        # return NaN, rand_n
        throw(ErrorException("rates_sum equal to zero in randchoose()"))
    end
end

function rejectionSample(rates::AbstractVector{Float64}, r_max::Float64)
    indices = eachindex(rates)
    while true
        i = rand(indices)
        if rates[i] > rand() * r_max
            return i
        end
    end
end

# Model entity sampling:
function choosePathogen(host_idx::Int64, population_idx::Int64, weight::Int64, model::Model, rand_n::Float64)
    return randChoose(
        rand_n,
        @views(model.populations[population_idx].hosts[host_idx].pathogen_weights[weight, :]),
        model.populations[population_idx].host_weights[weight, host_idx],
        # assumes population's host weights are the sum of all pathogen weights, should be true
        regenerate_rand=true
    )
end

function choosePathogenFromPopulationFraction(host_idx::Int64, population_idx::Int64, model::Model, rand_n::Float64)
    return randChoose(
        rand_n,
        @views(model.populations[population_idx].hosts[host_idx].pathogen_fractions),
        1.0,
        regenerate_rand=true
    )
end

function chooseResponse(host_idx::Int64, population_idx::Int64, weight::Int64, model::Model, rand_n::Float64)
    return randChoose(
        rand_n,
        @views(model.populations[population_idx].hosts[host_idx].response_weights[weight-RESPONSE_EVENTS[1]+1, :]),
        model.populations[population_idx].host_weights[weight, host_idx],
        # assumes population's host weights are the sum of all response weights, should be true
        regenerate_rand=true
    )
end

function chooseHost(population_idx::Int64, weight::Int64, model::Model, rand_n::Float64)
    if weight > NUM_EVENTS
        return sample(model.populations[population_idx].host_weights_receive_with_coefficient_samplers[weight-CHOICE_MODIFIERS[1]+1]), rand_n
        # return randChoose(
        #     rand_n,
        #     @views(model.populations[population_idx].host_weights_receive_with_coefficient[weight-CHOICE_MODIFIERS[1]+1, :]),
        #     model.population_weights_receive[weight-CHOICE_MODIFIERS[1]+1, population_idx],
        #     #TODO: revise, maybe this is another field?
        #     regenerate_rand=true
        # )
    else
        return sample(model.populations[population_idx].host_weights_with_coefficient_samplers[weight]), rand_n
        # return randChoose(
        #     rand_n,
        #     @views(model.populations[population_idx].host_weights_with_coefficient[weight, :]),
        #     model.population_weights[weight, population_idx],
        #     #TODO: revise, maybe this is another field?
        #     regenerate_rand=true
        # )
    end
end

function choosePopulation(weight::Int64, model::Model, rand_n::Float64)
    if weight > NUM_EVENTS
        return randChoose(
            rand_n,
            @views(model.population_weights_receive[weight-CHOICE_MODIFIERS[1]+1, :]),
            model.population_weights_receive_sums[weight-CHOICE_MODIFIERS[1]+1],
            regenerate_rand=true
        )
    else
        return randChoose(
            rand_n,
            @views(model.population_weights[weight, :]),
            model.event_rates[weight],
            regenerate_rand=true
        )
    end
end

function choosePopulationReceiveContact(emitter_idx::Int64, model::Model, rand_n::Float64)
    return randChoose(
        rand_n,
        @views(model.population_contact_weights_receive[:, emitter_idx]),
        model.population_contact_weights_receive_sums[emitter_idx],
        regenerate_rand=true
    )
end

function choosePopulationReceiveTransition(emitter_idx::Int64, model::Model, rand_n::Float64)
    return randChoose(
        rand_n,
        @views(model.population_transition_weights_receive[:, emitter_idx]),
        model.population_transition_weights_receive_sums[emitter_idx],
        regenerate_rand=true
    )
end

function chooseEvent(model::Model, rand_n::Float64)
    return randChoose(
        rand_n,
        model.event_rates,
        model.event_rates_sum,
        regenerate_rand=true
    )
end

function choosePathogen(weight::Int64, model::Model, rand_n::Float64)
    pop_idx, rand_n = choosePopulation(weight, model, rand_n)
    host_idx, rand_n = chooseHost(pop_idx, weight, model, rand_n)
    pathogen_idx, rand_n = choosePathogen(host_idx, pop_idx, weight, model, rand_n)

    return pathogen_idx, host_idx, pop_idx, rand_n
end

function chooseResponse(weight::Int64, model::Model, rand_n::Float64)
    pop_idx, rand_n = choosePopulation(weight, model, rand_n)
    host_idx, rand_n = chooseHost(pop_idx, weight, model, rand_n)
    response_idx, rand_n = chooseResponse(host_idx, pop_idx, weight, model, rand_n)

    return response_idx, host_idx, pop_idx, rand_n
end

function chooseHost(weight::Int64, model::Model, rand_n::Float64)
    pop_idx, rand_n = choosePopulation(weight, model, rand_n)
    host_idx, rand_n = chooseHost(pop_idx, weight, model, rand_n)

    return host_idx, pop_idx, rand_n
end
