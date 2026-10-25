using StaticArrays
using Flexle

# Model entity initializers

function newPathogen!(
    sequence::String, population::Population, type::PathogenType;
    parents::MVector{2,Union{Pathogen,Nothing}}=MVector{2,Union{Pathogen,Nothing}}([nothing, nothing]),
    birth_time::Float64=0.0)
    population.pathogens[sequence] = Pathogen(
        parents, birth_time, sequence,
        pathogenSequenceSpecificCoefficients(sequence, population.id, type),
        pathogenSequenceHostwideCoefficients(sequence, population.id, type),
        type
    )

    return population.pathogens[sequence]
end

function newResponse!(
    imprinted_pathogen::Union{Pathogen,Nothing}, matured_pathogen::Union{Pathogen,Nothing}, host_sequence::String,
    existing_responses::Dict{Tuple{String,String,String,String},Response}, population_id::String, type::ResponseType;
    parents::MVector{2,Union{Response,Nothing}}=MVector{2,Union{Response,Nothing}}([nothing, nothing]),
    birth_time::Float64=0.0)
    isnothing(imprinted_pathogen) ? imprinted_seq = "" : imprinted_seq = imprinted_pathogen.sequence
    isnothing(matured_pathogen) ? matured_seq = "" : matured_seq = matured_pathogen.sequence

    existing_responses[(imprinted_seq, matured_seq, host_sequence, type.id)] = Response(
        parents, birth_time, host_sequence, imprinted_pathogen, matured_pathogen,
        responseStaticSpecificCoefficients(imprinted_seq, matured_seq, host_sequence, population_id, type),
        responseStaticHostwideCoefficients(imprinted_seq, matured_seq, host_sequence, population_id, type),
        type
    )

    return existing_responses[(imprinted_seq, matured_seq, host_sequence, type.id)]
end

function newHost!(sequence::String, type::HostType, population::Population, model::Model;
    parents::MVector{2,Union{Host,Nothing}}=MVector{2,Union{Host,Nothing}}([nothing, nothing]),
    birth_time::Float64=0.0)
    addHostToPopulation!(
        Host(
            length(population.hosts) + 1,
            parents,
            birth_time,
            sequence,
            Vector{Pathogen}(undef, 0), Vector{Response}(undef, 0),
            Vector{Float64}(undef, 0),
            Dict{String,MVector{3,Float64}}(),
            Matrix{Float64}(undef, NUM_PATHOGEN_EVENTS, 0),
            Matrix{Float64}(undef, NUM_RESPONSE_EVENTS, 0),
            hostSequenceCoefficients(sequence, population.id, type),
            type,
        ),
        population, model
    )

    return population.hosts[end]
end

function staticHost(host::Host)
    return StaticHost(host.id, host.sequence, copy(host.pathogens), copy(host.responses))
end

function newPopulation!(id::String, parameters::PopulationType, model::Model)
    push!(model.populations, Population(
        id, parameters,
        Dict{String,Pathogen}(), Dict{Tuple{String,String,String,String},Response}(),
        Vector{Host}(undef, 0),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS, 0),
        # Matrix{Float64}(undef, NUM_EVENTS, 0),
        # Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS, 0),
        # MVector{NUM_EVENTS,FlexleSampler}(repeat([FlexleSampler()], outer=NUM_EVENTS)),
        MVector{NUM_EVENTS,FlexleSampler}(flexleSamplers(NUM_EVENTS)),
        # MVector{NUM_CHOICE_MODIFIERS,FlexleSampler}(repeat([FlexleSampler()], outer=NUM_CHOICE_MODIFIERS)),
        MVector{NUM_CHOICE_MODIFIERS,FlexleSampler}(flexleSamplers(NUM_CHOICE_MODIFIERS)),
        0.0, 0.0,
        zeros(Float64, length(model.populations)),
        zeros(Float64, length(model.populations)),
        zeros(Int64, NUM_COMPARTMENTS),
        zeros(Int64, NUM_SAMPLING_COEFFICIENTS-1),
    ))
    model.population_dict[id] = length(model.populations)
    model.population_weights = catCol(model.population_weights, zeros(Float64, NUM_EVENTS))
    model.population_weights_receive = catCol(model.population_weights_receive, zeros(Float64, NUM_CHOICE_MODIFIERS))
    model.population_weights_receive[ RECEIVE_TRANSITION-CHOICE_MODIFIERS[1]+1, model.population_dict[id] ] = parameters.base_coefficients[RECEIVE_TRANSITION]

    for pop in model.populations
        push!(pop.population_contact_coefficients, 0.0)
        push!(pop.population_transition_coefficients, 0.0)
    end

    if length(model.populations) > 1 # if more than just this population we just added
        new_weights = zeros( Float64, length(model.populations), length(model.populations) )
        new_weights[1:length(model.populations)-1, 1:length(model.populations)-1] = model.population_contact_weights_receive
        model.population_contact_weights_receive = copy(new_weights)

        new_weights[1:length(model.populations)-1, 1:length(model.populations)-1] = model.population_transition_weights_receive
        model.population_transition_weights_receive = copy(new_weights)
    else
        model.population_contact_weights_receive = Matrix{Float64}(undef, 1, 1)
        model.population_transition_weights_receive = Matrix{Float64}(undef, 1, 1)
    end

    push!(model.population_contact_weights_receive_sums, 0.0)
    push!(model.population_transition_weights_receive_sums, 0.0)

    # By default, self contact is set to 1.0
    setPopulationContactCoefficient!(
        model.populations[length(model.populations)], model.populations[length(model.populations)], 1.0, model
    )

    return model.populations[end]
end

function newModel()
    return Model(
        Vector{Population}(undef, 0),
        Dict{String,Int64}(),
        Matrix{Float64}(undef, NUM_EVENTS, 0),
        Matrix{Float64}(undef, NUM_CHOICE_MODIFIERS, 0),
        zeros(SVector{NUM_CHOICE_MODIFIERS,Float64}),
        Matrix{Float64}(undef, 0, 0),
        Matrix{Float64}(undef, 0, 0),
        Vector{Float64}(undef, 0),
        Vector{Float64}(undef, 0),
        zeros(SVector{NUM_EVENTS,Float64}),
        0.0,
        zeros(Int64, NUM_SAMPLING_COEFFICIENTS-1),
        0.0,
    )
end
