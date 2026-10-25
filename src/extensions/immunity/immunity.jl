function weightedInteractionPathogenArithmeticMean(pathogen::Pathogen, host::Host, population_id::String, evt::Int64)
    # reactivity-weighted arithmetic mean of interaction coefficients
    if length(host.responses) > 0
        reac_sum = 0.0
        numerator_sum = 0.0
        reac = 0.0
        for response in host.responses
            reac = reactivityCoefficient(pathogen, response, host, population_id)
            reac_sum += reac
            numerator_sum += reac * responseInteractionSpecificCoefficient(pathogen, response, host, population_id, evt)
        end

        return numerator_sum / reac_sum
    else
        return 1.0
    end
end

# function transmissionEfficiencyArithmeticMean(pathogen::Pathogen, host::Host)
#     # reactivity-weighted arithmetic mean of infection coefficients
#     if length(host.responses) > 0
#         reac_sum = 0.0
#         numerator_sum = 0.0
#         reac = 0.0
#         for response in host.responses
#             reac = reactivityCoefficient(pathogen, response, host)
#             reac_sum += reac
#             numerator_sum += reac * responseInteractionSpecificCoefficient(
#                 pathogen, response, host, TRANSMISSION_EFFICIENCY
#             )
#         end

#         return numerator_sum / reac_sum
#     else
#         return 1.0
#     end
# end

function weightedInteractionPathogenWinnerTakesAll(pathogen::Pathogen, host::Host, population_id::String, evt::Int64)
    # In case of ties, takes response that was acquired first
    if length(host.responses) > 0
        dominant_reaction = 1.0
        dominant_reactivity = 0.0
        for response in host.responses
            reactivity = reactivityCoefficient(pathogen, response, host, population_id)
            if reactivity > dominant_reactivity
                dominant_reaction = responseInteractionSpecificCoefficient(pathogen, response, host, population_id, evt)
                dominant_reactivity = reactivity
            end
        end

        return dominant_reaction
    else
        return 1.0
    end
end

function weightedInteractionResponseWinnerTakesAll(response::Response, host::Host, population_id::String, evt::Int64)
    # In case of ties, takes response that was acquired first
    if length(host.pathogens) > 0
        dominant_reaction = 1.0
        dominant_reactivity = 0.0
        for pathogen in host.pathogens
            reactivity = reactivityCoefficient(pathogen, response, host, population_id)
            if reactivity > dominant_reactivity
                dominant_reaction = responseInteractionSpecificCoefficient(pathogen, response, host, population_id, evt)
                dominant_reactivity = reactivity
            end
        end

        return dominant_reaction
    else
        return 1.0
    end
end

function weightedInteractionHostwideProduct(host::Host, population_id::String, evt::Int64)
    reaction = 1.0
    if length(host.pathogens) > 0
        if length(host.responses) > 0
            for pathogen in host.pathogens
                for response in host.responses
                    reaction = reaction * responseInteractionHostwideCoefficient(pathogen, response, host, population_id, evt)
                end
            end
        end
    end

    return reaction
end

# function transmissionEfficiencyWinnerTakesAll(pathogen::Pathogen, host::Host)
#     if length(host.responses) > 0
#         dominant_reaction = 1.0
#         dominant_reactivity = 0.0
#         for response in host.responses
#             reactivity = reactivityCoefficient(pathogen, response, host)
#             if reactivity > dominant_reactivity
#                 dominant_reaction = reactivity * responseInteractionSpecificCoefficient(
#                     pathogen, response, host, TRANSMISSION_EFFICIENCY
#                 )
#             end
#         end

#         return dominant_reaction
#     else
#         return 1.0
#     end
# end

function deNovoResponse(
    pathogen::Pathogen, host::Host, population_id::String,
    existing_responses::Dict{Tuple{String,String,String,String},Response},
    response_types::Dict{String,ResponseType}, birth_time::Float64; response_type_id::String="Default")
    if haskey(existing_responses, (host.sequence, pathogen.sequence, "", response_type_id))
        return existing_responses[(host.sequence, pathogen.sequence, "", response_type_id)]
    else
        return newResponse!(
            pathogen, nothing, host.sequence, existing_responses, population_id, response_types[response_type_id],
            parents=MVector{2,Union{Response,Nothing}}([nothing, nothing]), birth_time=birth_time
        )
    end
end
