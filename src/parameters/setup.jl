using StaticArrays

# Model parameter setup initializers

function newPathogenType(
    id::String;
    template::PathogenType=DEFAULT_PATHOGEN_TYPE,
    num_loci::Union{Nothing,Int64}=nothing,
    possible_alleles::Union{Nothing,String}=nothing,

    # Each element takes seq argument, returns Float64
    pathogenEstablishmentSpecificCoefficient::Union{Nothing,Function}=nothing,
    clearanceSpecificCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionSpecificCoefficient::Union{Nothing,Function}=nothing,
    contactSpecificCoefficient::Union{Nothing,Function}=nothing,
    intrahostFitnessSpecificCoefficient::Union{Nothing,Function}=nothing,
    generationsPerTransmissionSpecificCoefficient::Union{Nothing,Function}=nothing,
    mutationsPerGenerationSpecificCoefficient::Union{Nothing,Function}=nothing,
    recombinationsPerGenerationSpecificCoefficient::Union{Nothing,Function}=nothing,
    inoculumSpecificCoefficient::Union{Nothing,Function}=nothing,
    transmissionEfficiencySpecificCoefficient::Union{Nothing,Function}=nothing,
    verticalTransmissionSpecificCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionUponClearanceSpecificCoefficient::Union{Nothing,Function}=nothing,

    # Each element takes seq argument, returns Float64
    pathogenEstablishmentHostwideCoefficient::Union{Nothing,Function}=nothing,
    clearanceHostwideCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionHostwideCoefficient::Union{Nothing,Function}=nothing,
    contactHostwideCoefficient::Union{Nothing,Function}=nothing,
    responseLossHostwideCoefficient::Union{Nothing,Function}=nothing,
    birthHostwideCoefficient::Union{Nothing,Function}=nothing,
    deathHostwideCoefficient::Union{Nothing,Function}=nothing,
    transitionHostwideCoefficient::Union{Nothing,Function}=nothing,
    receiveTransitionHostwideCoefficient::Union{Nothing,Function}=nothing,
    receiveContactHostwideCoefficient::Union{Nothing,Function}=nothing,
    intrahostFitnessHostwideCoefficient::Union{Nothing,Function}=nothing,
    generationsPerTransmissionHostwideCoefficient::Union{Nothing,Function}=nothing,
    mutationsPerGenerationHostwideCoefficient::Union{Nothing,Function}=nothing,
    recombinationsPerGenerationHostwideCoefficient::Union{Nothing,Function}=nothing,
    inoculumHostwideCoefficient::Union{Nothing,Function}=nothing,
    transmissionEfficiencyHostwideCoefficient::Union{Nothing,Function}=nothing,
    verticalTransmissionHostwideCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionUponClearanceHostwideCoefficient::Union{Nothing,Function}=nothing,
    responseInheritanceHostwideCoefficient::Union{Nothing,Function}=nothing,
    hostMutationsUponBirthHostwideCoefficient::Union{Nothing,Function}=nothing,
    hostRecombinationsUponBirthHostwideCoefficient::Union{Nothing,Function}=nothing,
)

    isnothing(num_loci) ? num_loci = template.num_loci : num_loci = num_loci
    isnothing(possible_alleles) ? possible_alleles = template.possible_alleles : possible_alleles = possible_alleles

    isnothing(pathogenEstablishmentSpecificCoefficient) ? pathogenEstablishmentSpecificCoefficient = template.specific_coefficient_functions[PATHOGEN_ESTABLISHMENT] : pathogenEstablishmentSpecificCoefficient = pathogenEstablishmentSpecificCoefficient
    isnothing(clearanceSpecificCoefficient) ? clearanceSpecificCoefficient = template.specific_coefficient_functions[CLEARANCE] : clearanceSpecificCoefficient = clearanceSpecificCoefficient
    isnothing(responseAcquisitionSpecificCoefficient) ? responseAcquisitionSpecificCoefficient = template.specific_coefficient_functions[RESPONSE_ACQUISITION] : responseAcquisitionSpecificCoefficient = responseAcquisitionSpecificCoefficient
    isnothing(contactSpecificCoefficient) ? contactSpecificCoefficient = template.specific_coefficient_functions[CONTACT] : contactSpecificCoefficient = contactSpecificCoefficient
    isnothing(intrahostFitnessSpecificCoefficient) ? intrahostFitnessSpecificCoefficient = template.specific_coefficient_functions[INTRAHOST_FITNESS] : intrahostFitnessSpecificCoefficient = intrahostFitnessSpecificCoefficient
    isnothing(generationsPerTransmissionSpecificCoefficient) ? generationsPerTransmissionSpecificCoefficient = template.specific_coefficient_functions[GENERATIONS_PER_TRANSMISSION] : generationsPerTransmissionSpecificCoefficient = generationsPerTransmissionSpecificCoefficient
    isnothing(mutationsPerGenerationSpecificCoefficient) ? mutationsPerGenerationSpecificCoefficient = template.specific_coefficient_functions[MUTATIONS_PER_GENERATION] : mutationsPerGenerationSpecificCoefficient = mutationsPerGenerationSpecificCoefficient
    isnothing(recombinationsPerGenerationSpecificCoefficient) ? recombinationsPerGenerationSpecificCoefficient = template.specific_coefficient_functions[RECOMBINATIONS_PER_GENERATION] : recombinationsPerGenerationSpecificCoefficient = recombinationsPerGenerationSpecificCoefficient
    isnothing(inoculumSpecificCoefficient) ? inoculumSpecificCoefficient = template.specific_coefficient_functions[INOCULUM] : inoculumSpecificCoefficient = inoculumSpecificCoefficient
    isnothing(transmissionEfficiencySpecificCoefficient) ? transmissionEfficiencySpecificCoefficient = template.specific_coefficient_functions[TRANSMISSION_EFFICIENCY] : transmissionEfficiencySpecificCoefficient = transmissionEfficiencySpecificCoefficient
    isnothing(verticalTransmissionSpecificCoefficient) ? verticalTransmissionSpecificCoefficient = template.specific_coefficient_functions[VERTICAL_TRANSMISSION] : verticalTransmissionSpecificCoefficient = verticalTransmissionSpecificCoefficient
    isnothing(responseAcquisitionUponClearanceSpecificCoefficient) ? responseAcquisitionUponClearanceSpecificCoefficient = template.specific_coefficient_functions[RESPONSE_ACQUISITION_UPON_CLEARANCE] : responseAcquisitionUponClearanceSpecificCoefficient = responseAcquisitionUponClearanceSpecificCoefficient

    isnothing(pathogenEstablishmentHostwideCoefficient) ? pathogenEstablishmentHostwideCoefficient = template.hostwide_coefficient_functions[PATHOGEN_ESTABLISHMENT] : pathogenEstablishmentHostwideCoefficient = pathogenEstablishmentHostwideCoefficient
    isnothing(clearanceHostwideCoefficient) ? clearanceHostwideCoefficient = template.hostwide_coefficient_functions[CLEARANCE] : clearanceHostwideCoefficient = clearanceHostwideCoefficient
    isnothing(responseAcquisitionHostwideCoefficient) ? responseAcquisitionHostwideCoefficient = template.hostwide_coefficient_functions[RESPONSE_ACQUISITION] : responseAcquisitionHostwideCoefficient = responseAcquisitionHostwideCoefficient
    isnothing(contactHostwideCoefficient) ? contactHostwideCoefficient = template.hostwide_coefficient_functions[CONTACT] : contactHostwideCoefficient = contactHostwideCoefficient
    isnothing(responseLossHostwideCoefficient) ? responseLossHostwideCoefficient = template.hostwide_coefficient_functions[RESPONSE_LOSS] : responseLossHostwideCoefficient = responseLossHostwideCoefficient
    isnothing(birthHostwideCoefficient) ? birthHostwideCoefficient = template.hostwide_coefficient_functions[BIRTH] : birthHostwideCoefficient = birthHostwideCoefficient
    isnothing(deathHostwideCoefficient) ? deathHostwideCoefficient = template.hostwide_coefficient_functions[DEATH] : deathHostwideCoefficient = deathHostwideCoefficient
    isnothing(transitionHostwideCoefficient) ? transitionHostwideCoefficient = template.hostwide_coefficient_functions[TRANSITION] : transitionHostwideCoefficient = transitionHostwideCoefficient
    isnothing(receiveTransitionHostwideCoefficient) ? receiveTransitionHostwideCoefficient = template.hostwide_coefficient_functions[RECEIVE_TRANSITION] : receiveTransitionHostwideCoefficient = receiveTransitionHostwideCoefficient
    isnothing(receiveContactHostwideCoefficient) ? receiveContactHostwideCoefficient = template.hostwide_coefficient_functions[RECEIVE_CONTACT] : receiveContactHostwideCoefficient = receiveContactHostwideCoefficient
    isnothing(intrahostFitnessHostwideCoefficient) ? intrahostFitnessHostwideCoefficient = template.hostwide_coefficient_functions[INTRAHOST_FITNESS] : intrahostFitnessHostwideCoefficient = intrahostFitnessHostwideCoefficient
    isnothing(generationsPerTransmissionHostwideCoefficient) ? generationsPerTransmissionHostwideCoefficient = template.hostwide_coefficient_functions[GENERATIONS_PER_TRANSMISSION] : generationsPerTransmissionHostwideCoefficient = generationsPerTransmissionHostwideCoefficient
    isnothing(mutationsPerGenerationHostwideCoefficient) ? mutationsPerGenerationHostwideCoefficient = template.hostwide_coefficient_functions[MUTATIONS_PER_GENERATION] : mutationsPerGenerationHostwideCoefficient = mutationsPerGenerationHostwideCoefficient
    isnothing(recombinationsPerGenerationHostwideCoefficient) ? recombinationsPerGenerationHostwideCoefficient = template.hostwide_coefficient_functions[RECOMBINATIONS_PER_GENERATION] : recombinationsPerGenerationHostwideCoefficient = recombinationsPerGenerationHostwideCoefficient
    isnothing(inoculumHostwideCoefficient) ? inoculumHostwideCoefficient = template.hostwide_coefficient_functions[INOCULUM] : inoculumHostwideCoefficient = inoculumHostwideCoefficient
    isnothing(transmissionEfficiencyHostwideCoefficient) ? transmissionEfficiencyHostwideCoefficient = template.hostwide_coefficient_functions[TRANSMISSION_EFFICIENCY] : transmissionEfficiencyHostwideCoefficient = transmissionEfficiencyHostwideCoefficient
    isnothing(verticalTransmissionHostwideCoefficient) ? verticalTransmissionHostwideCoefficient = template.hostwide_coefficient_functions[VERTICAL_TRANSMISSION] : verticalTransmissionHostwideCoefficient = verticalTransmissionHostwideCoefficient
    isnothing(responseAcquisitionUponClearanceHostwideCoefficient) ? responseAcquisitionUponClearanceHostwideCoefficient = template.hostwide_coefficient_functions[RESPONSE_ACQUISITION_UPON_CLEARANCE] : responseAcquisitionUponClearanceHostwideCoefficient = responseAcquisitionUponClearanceHostwideCoefficient
    isnothing(responseInheritanceHostwideCoefficient) ? responseInheritanceHostwideCoefficient = template.hostwide_coefficient_functions[RESPONSE_INHERITANCE] : responseInheritanceHostwideCoefficient = responseInheritanceHostwideCoefficient
    isnothing(hostMutationsUponBirthHostwideCoefficient) ? hostMutationsUponBirthHostwideCoefficient = template.hostwide_coefficient_functions[HOST_MUTATIONS_UPON_BIRTH] : hostMutationsUponBirthHostwideCoefficient = hostMutationsUponBirthHostwideCoefficient
    isnothing(hostRecombinationsUponBirthHostwideCoefficient) ? hostRecombinationsUponBirthHostwideCoefficient = template.hostwide_coefficient_functions[HOST_RECOMBINATIONS_UPON_BIRTH] : hostRecombinationsUponBirthHostwideCoefficient = hostRecombinationsUponBirthHostwideCoefficient


    return PathogenType(
        id,
        num_loci,
        possible_alleles,
        SA[ # order defined in COEFFICIENTS
            pathogenEstablishmentSpecificCoefficient, clearanceSpecificCoefficient, responseAcquisitionSpecificCoefficient, contactSpecificCoefficient,
            (g::String,pop_id::String)->1.0, (g::String,pop_id::String)->1.0, (g::String,pop_id::String)->1.0,
            (g::String,pop_id::String)->1.0, (g::String,pop_id::String)->1.0, (g::String,pop_id::String)->1.0,
            intrahostFitnessSpecificCoefficient, generationsPerTransmissionSpecificCoefficient,
            mutationsPerGenerationSpecificCoefficient, recombinationsPerGenerationSpecificCoefficient, inoculumSpecificCoefficient,
            transmissionEfficiencySpecificCoefficient, verticalTransmissionSpecificCoefficient, responseAcquisitionUponClearanceSpecificCoefficient,
            (g::String,pop_id::String)->1.0, (g::String,pop_id::String)->1.0, (g::String,pop_id::String)->1.0,
        ],
        SA[ # order defined in COEFFICIENTS
            pathogenEstablishmentHostwideCoefficient, clearanceHostwideCoefficient, responseAcquisitionHostwideCoefficient,
            contactHostwideCoefficient, responseLossHostwideCoefficient,
            birthHostwideCoefficient, deathHostwideCoefficient, transitionHostwideCoefficient,
            receiveTransitionHostwideCoefficient, receiveContactHostwideCoefficient, intrahostFitnessHostwideCoefficient,
            generationsPerTransmissionHostwideCoefficient,
            mutationsPerGenerationHostwideCoefficient, recombinationsPerGenerationHostwideCoefficient, inoculumHostwideCoefficient,
            transmissionEfficiencyHostwideCoefficient, verticalTransmissionHostwideCoefficient, responseAcquisitionUponClearanceHostwideCoefficient,
            responseInheritanceHostwideCoefficient, hostMutationsUponBirthHostwideCoefficient, hostRecombinationsUponBirthHostwideCoefficient,
        ],
    )
end

function newResponseType(
    id::String;
    template::ResponseType=DEFAULT_RESPONSE_TYPE,
    reactivityCoefficient::Union{Nothing,Function}=nothing,
    # takes host, imprinted, matured, and infecting sequences and returns Float64 coefficient

    # Each takes host, imprinted, matured sequences and returns Float64 coefficient
    responseLossStaticSpecificCoefficient::Union{Nothing,Function}=nothing,
    responseInheritanceStaticSpecificCoefficient::Union{Nothing,Function}=nothing,

    # Each takes host, imprinted, matured, and infecting sequences and returns Float64 coefficient
    pathogenEstablishmentInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    clearanceInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    contactInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    responseLossInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    birthInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    deathInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    transitionInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    receiveTransitionInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    receiveContactInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    intrahostFitnessInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    generationsPerTransmissionInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    mutationsPerGenerationInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    recombinationsPerGenerationInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    inoculumInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    transmissionEfficiencyInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    verticalTransmissionInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionUponClearanceInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,
    responseInheritanceInteractionSpecificCoefficient::Union{Nothing,Function}=nothing,

    # Each takes host, imprinted, matured sequences and returns Float64 coefficient
    pathogenEstablishmentStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    clearanceStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    contactStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    responseLossStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    birthStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    deathStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    transitionStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    receiveTransitionStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    receiveContactStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    intrahostFitnessStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    generationsPerTransmissionStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    mutationsPerGenerationStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    recombinationsPerGenerationStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    inoculumStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    transmissionEfficiencyStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    verticalTransmissionStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionUponClearanceStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    responseInheritanceStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    hostMutationsUponBirthStaticHostwideCoefficient::Union{Nothing,Function}=nothing,
    hostRecombinationsUponBirthStaticHostwideCoefficient::Union{Nothing,Function}=nothing,

    # Each takes host, imprinted, matured, and infecting sequences and returns Float64 coefficient
    pathogenEstablishmentInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    clearanceInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    contactInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    responseLossInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    birthInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    deathInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    transitionInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    receiveTransitionInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    receiveContactInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    intrahostFitnessInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    generationsPerTransmissionInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    mutationsPerGenerationInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    recombinationsPerGenerationInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    inoculumInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    transmissionEfficiencyInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    verticalTransmissionInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionUponClearanceInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    responseInheritanceInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    hostMutationsUponBirthInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
    hostRecombinationsUponBirthInteractionHostwideCoefficient::Union{Nothing,Function}=nothing,
)

    isnothing(reactivityCoefficient) ? reactivityCoefficient = template.reactivityCoefficient : reactivityCoefficient = reactivityCoefficient

    isnothing(responseLossStaticSpecificCoefficient) ? responseLossStaticSpecificCoefficient = template.static_specific_coefficient_functions[RESPONSE_LOSS] : responseLossStaticSpecificCoefficient = responseLossStaticSpecificCoefficient
    isnothing(responseInheritanceStaticSpecificCoefficient) ? responseInheritanceStaticSpecificCoefficient = template.static_specific_coefficient_functions[RESPONSE_INHERITANCE] : responseInheritanceStaticSpecificCoefficient = responseInheritanceStaticSpecificCoefficient

    isnothing(pathogenEstablishmentInteractionSpecificCoefficient) ? pathogenEstablishmentInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[PATHOGEN_ESTABLISHMENT] : pathogenEstablishmentInteractionSpecificCoefficient = pathogenEstablishmentInteractionSpecificCoefficient
    isnothing(clearanceInteractionSpecificCoefficient) ? clearanceInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[CLEARANCE] : clearanceInteractionSpecificCoefficient = clearanceInteractionSpecificCoefficient
    isnothing(responseAcquisitionInteractionSpecificCoefficient) ? responseAcquisitionInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[RESPONSE_ACQUISITION] : responseAcquisitionInteractionSpecificCoefficient = responseAcquisitionInteractionSpecificCoefficient
    isnothing(contactInteractionSpecificCoefficient) ? contactInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[CONTACT] : contactInteractionSpecificCoefficient = contactInteractionSpecificCoefficient
    isnothing(responseLossInteractionSpecificCoefficient) ? responseLossInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[RESPONSE_LOSS] : responseLossInteractionSpecificCoefficient = responseLossInteractionSpecificCoefficient
    isnothing(birthInteractionSpecificCoefficient) ? birthInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[BIRTH] : birthInteractionSpecificCoefficient = birthInteractionSpecificCoefficient
    isnothing(deathInteractionSpecificCoefficient) ? deathInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[DEATH] : deathInteractionSpecificCoefficient = deathInteractionSpecificCoefficient
    isnothing(transitionInteractionSpecificCoefficient) ? transitionInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[TRANSITION] : transitionInteractionSpecificCoefficient = transitionInteractionSpecificCoefficient
    isnothing(receiveTransitionInteractionSpecificCoefficient) ? receiveTransitionInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[RECEIVE_TRANSITION] : receiveTransitionInteractionSpecificCoefficient = receiveTransitionInteractionSpecificCoefficient
    isnothing(receiveContactInteractionSpecificCoefficient) ? receiveContactInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[RECEIVE_CONTACT] : receiveContactInteractionSpecificCoefficient = receiveContactInteractionSpecificCoefficient
    isnothing(intrahostFitnessInteractionSpecificCoefficient) ? intrahostFitnessInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[INTRAHOST_FITNESS] : intrahostFitnessInteractionSpecificCoefficient = intrahostFitnessInteractionSpecificCoefficient
    isnothing(generationsPerTransmissionInteractionSpecificCoefficient) ? generationsPerTransmissionInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[GENERATIONS_PER_TRANSMISSION] : generationsPerTransmissionInteractionSpecificCoefficient = generationsPerTransmissionInteractionSpecificCoefficient
    isnothing(mutationsPerGenerationInteractionSpecificCoefficient) ? mutationsPerGenerationInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[MUTATIONS_PER_GENERATION] : mutationsPerGenerationInteractionSpecificCoefficient = mutationsPerGenerationInteractionSpecificCoefficient
    isnothing(recombinationsPerGenerationInteractionSpecificCoefficient) ? recombinationsPerGenerationInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[RECOMBINATIONS_PER_GENERATION] : recombinationsPerGenerationInteractionSpecificCoefficient = recombinationsPerGenerationInteractionSpecificCoefficient
    isnothing(inoculumInteractionSpecificCoefficient) ? inoculumInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[INOCULUM] : inoculumInteractionSpecificCoefficient = inoculumInteractionSpecificCoefficient
    isnothing(transmissionEfficiencyInteractionSpecificCoefficient) ? transmissionEfficiencyInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[TRANSMISSION_EFFICIENCY] : transmissionEfficiencyInteractionSpecificCoefficient = transmissionEfficiencyInteractionSpecificCoefficient
    isnothing(verticalTransmissionInteractionSpecificCoefficient) ? verticalTransmissionInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[VERTICAL_TRANSMISSION] : verticalTransmissionInteractionSpecificCoefficient = verticalTransmissionInteractionSpecificCoefficient
    isnothing(responseAcquisitionUponClearanceInteractionSpecificCoefficient) ? responseAcquisitionUponClearanceInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[RESPONSE_ACQUISITION_UPON_CLEARANCE] : responseAcquisitionUponClearanceInteractionSpecificCoefficient = responseAcquisitionUponClearanceInteractionSpecificCoefficient
    isnothing(responseInheritanceInteractionSpecificCoefficient) ? responseInheritanceInteractionSpecificCoefficient = template.interaction_specific_coefficient_functions[RESPONSE_INHERITANCE] : responseInheritanceInteractionSpecificCoefficient = responseInheritanceInteractionSpecificCoefficient

    isnothing(pathogenEstablishmentStaticHostwideCoefficient) ? pathogenEstablishmentStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[PATHOGEN_ESTABLISHMENT] : pathogenEstablishmentStaticHostwideCoefficient = pathogenEstablishmentStaticHostwideCoefficient
    isnothing(clearanceStaticHostwideCoefficient) ? clearanceStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[CLEARANCE] : clearanceStaticHostwideCoefficient = clearanceStaticHostwideCoefficient
    isnothing(responseAcquisitionStaticHostwideCoefficient) ? responseAcquisitionStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[RESPONSE_ACQUISITION] : responseAcquisitionStaticHostwideCoefficient = responseAcquisitionStaticHostwideCoefficient
    isnothing(contactStaticHostwideCoefficient) ? contactStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[CONTACT] : contactStaticHostwideCoefficient = contactStaticHostwideCoefficient
    isnothing(responseLossStaticHostwideCoefficient) ? responseLossStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[RESPONSE_LOSS] : responseLossStaticHostwideCoefficient = responseLossStaticHostwideCoefficient
    isnothing(birthStaticHostwideCoefficient) ? birthStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[BIRTH] : birthStaticHostwideCoefficient = birthStaticHostwideCoefficient
    isnothing(deathStaticHostwideCoefficient) ? deathStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[DEATH] : deathStaticHostwideCoefficient = deathStaticHostwideCoefficient
    isnothing(transitionStaticHostwideCoefficient) ? transitionStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[TRANSITION] : transitionStaticHostwideCoefficient = transitionStaticHostwideCoefficient
    isnothing(receiveTransitionStaticHostwideCoefficient) ? receiveTransitionStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[RECEIVE_TRANSITION] : receiveTransitionStaticHostwideCoefficient = receiveTransitionStaticHostwideCoefficient
    isnothing(receiveContactStaticHostwideCoefficient) ? receiveContactStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[RECEIVE_CONTACT] : receiveContactStaticHostwideCoefficient = receiveContactStaticHostwideCoefficient
    isnothing(intrahostFitnessStaticHostwideCoefficient) ? intrahostFitnessStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[INTRAHOST_FITNESS] : intrahostFitnessStaticHostwideCoefficient = intrahostFitnessStaticHostwideCoefficient
    isnothing(generationsPerTransmissionStaticHostwideCoefficient) ? generationsPerTransmissionStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[GENERATIONS_PER_TRANSMISSION] : generationsPerTransmissionStaticHostwideCoefficient = generationsPerTransmissionStaticHostwideCoefficient
    isnothing(mutationsPerGenerationStaticHostwideCoefficient) ? mutationsPerGenerationStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[MUTATIONS_PER_GENERATION] : mutationsPerGenerationStaticHostwideCoefficient = mutationsPerGenerationStaticHostwideCoefficient
    isnothing(recombinationsPerGenerationStaticHostwideCoefficient) ? recombinationsPerGenerationStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[RECOMBINATIONS_PER_GENERATION] : recombinationsPerGenerationStaticHostwideCoefficient = recombinationsPerGenerationStaticHostwideCoefficient
    isnothing(inoculumStaticHostwideCoefficient) ? inoculumStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[INOCULUM] : inoculumStaticHostwideCoefficient = inoculumStaticHostwideCoefficient
    isnothing(transmissionEfficiencyStaticHostwideCoefficient) ? transmissionEfficiencyStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[TRANSMISSION_EFFICIENCY] : transmissionEfficiencyStaticHostwideCoefficient = transmissionEfficiencyStaticHostwideCoefficient
    isnothing(verticalTransmissionStaticHostwideCoefficient) ? verticalTransmissionStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[VERTICAL_TRANSMISSION] : verticalTransmissionStaticHostwideCoefficient = verticalTransmissionStaticHostwideCoefficient
    isnothing(responseAcquisitionUponClearanceStaticHostwideCoefficient) ? responseAcquisitionUponClearanceStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[RESPONSE_ACQUISITION_UPON_CLEARANCE] : responseAcquisitionUponClearanceStaticHostwideCoefficient = responseAcquisitionUponClearanceStaticHostwideCoefficient
    isnothing(responseInheritanceStaticHostwideCoefficient) ? responseInheritanceStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[RESPONSE_INHERITANCE] : responseInheritanceStaticHostwideCoefficient = responseInheritanceStaticHostwideCoefficient
    isnothing(hostMutationsUponBirthStaticHostwideCoefficient) ? hostMutationsUponBirthStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[HOST_MUTATIONS_UPON_BIRTH] : hostMutationsUponBirthStaticHostwideCoefficient = hostMutationsUponBirthStaticHostwideCoefficient
    isnothing(hostRecombinationsUponBirthStaticHostwideCoefficient) ? hostRecombinationsUponBirthStaticHostwideCoefficient = template.static_hostwide_coefficient_functions[HOST_RECOMBINATIONS_UPON_BIRTH] : hostRecombinationsUponBirthStaticHostwideCoefficient = hostRecombinationsUponBirthStaticHostwideCoefficient

    isnothing(pathogenEstablishmentInteractionHostwideCoefficient) ? pathogenEstablishmentInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[PATHOGEN_ESTABLISHMENT] : pathogenEstablishmentInteractionHostwideCoefficient = pathogenEstablishmentInteractionHostwideCoefficient
    isnothing(clearanceInteractionHostwideCoefficient) ? clearanceInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[CLEARANCE] : clearanceInteractionHostwideCoefficient = clearanceInteractionHostwideCoefficient
    isnothing(responseAcquisitionInteractionHostwideCoefficient) ? responseAcquisitionInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[RESPONSE_ACQUISITION] : responseAcquisitionInteractionHostwideCoefficient = responseAcquisitionInteractionHostwideCoefficient
    isnothing(contactInteractionHostwideCoefficient) ? contactInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[CONTACT] : contactInteractionHostwideCoefficient = contactInteractionHostwideCoefficient
    isnothing(responseLossInteractionHostwideCoefficient) ? responseLossInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[RESPONSE_LOSS] : responseLossInteractionHostwideCoefficient = responseLossInteractionHostwideCoefficient
    isnothing(birthInteractionHostwideCoefficient) ? birthInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[BIRTH] : birthInteractionHostwideCoefficient = birthInteractionHostwideCoefficient
    isnothing(deathInteractionHostwideCoefficient) ? deathInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[DEATH] : deathInteractionHostwideCoefficient = deathInteractionHostwideCoefficient
    isnothing(transitionInteractionHostwideCoefficient) ? transitionInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[TRANSITION] : transitionInteractionHostwideCoefficient = transitionInteractionHostwideCoefficient
    isnothing(receiveTransitionInteractionHostwideCoefficient) ? receiveTransitionInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[RECEIVE_TRANSITION] : receiveTransitionInteractionHostwideCoefficient = receiveTransitionInteractionHostwideCoefficient
    isnothing(receiveContactInteractionHostwideCoefficient) ? receiveContactInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[RECEIVE_CONTACT] : receiveContactInteractionHostwideCoefficient = receiveContactInteractionHostwideCoefficient
    isnothing(intrahostFitnessInteractionHostwideCoefficient) ? intrahostFitnessInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[INTRAHOST_FITNESS] : intrahostFitnessInteractionHostwideCoefficient = intrahostFitnessInteractionHostwideCoefficient
    isnothing(generationsPerTransmissionInteractionHostwideCoefficient) ? generationsPerTransmissionInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[GENERATIONS_PER_TRANSMISSION] : generationsPerTransmissionInteractionHostwideCoefficient = generationsPerTransmissionInteractionHostwideCoefficient
    isnothing(mutationsPerGenerationInteractionHostwideCoefficient) ? mutationsPerGenerationInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[MUTATIONS_PER_GENERATION] : mutationsPerGenerationInteractionHostwideCoefficient = mutationsPerGenerationInteractionHostwideCoefficient
    isnothing(recombinationsPerGenerationInteractionHostwideCoefficient) ? recombinationsPerGenerationInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[RECOMBINATIONS_PER_GENERATION] : recombinationsPerGenerationInteractionHostwideCoefficient = recombinationsPerGenerationInteractionHostwideCoefficient
    isnothing(inoculumInteractionHostwideCoefficient) ? inoculumInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[INOCULUM] : inoculumInteractionHostwideCoefficient = inoculumInteractionHostwideCoefficient
    isnothing(transmissionEfficiencyInteractionHostwideCoefficient) ? transmissionEfficiencyInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[TRANSMISSION_EFFICIENCY] : transmissionEfficiencyInteractionHostwideCoefficient = transmissionEfficiencyInteractionHostwideCoefficient
    isnothing(verticalTransmissionInteractionHostwideCoefficient) ? verticalTransmissionInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[VERTICAL_TRANSMISSION] : verticalTransmissionInteractionHostwideCoefficient = verticalTransmissionInteractionHostwideCoefficient
    isnothing(responseAcquisitionUponClearanceInteractionHostwideCoefficient) ? responseAcquisitionUponClearanceInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[RESPONSE_ACQUISITION_UPON_CLEARANCE] : responseAcquisitionUponClearanceInteractionHostwideCoefficient = responseAcquisitionUponClearanceInteractionHostwideCoefficient
    isnothing(responseInheritanceInteractionHostwideCoefficient) ? responseInheritanceInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[RESPONSE_INHERITANCE] : responseInheritanceInteractionHostwideCoefficient = responseInheritanceInteractionHostwideCoefficient
    isnothing(hostMutationsUponBirthInteractionHostwideCoefficient) ? hostMutationsUponBirthInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[HOST_MUTATIONS_UPON_BIRTH] : hostMutationsUponBirthInteractionHostwideCoefficient = hostMutationsUponBirthInteractionHostwideCoefficient
    isnothing(hostRecombinationsUponBirthInteractionHostwideCoefficient) ? hostRecombinationsUponBirthInteractionHostwideCoefficient = template.interaction_hostwide_coefficient_functions[HOST_RECOMBINATIONS_UPON_BIRTH] : hostRecombinationsUponBirthInteractionHostwideCoefficient = hostRecombinationsUponBirthInteractionHostwideCoefficient

    return ResponseType(
        id,
        reactivityCoefficient,
        SA[ # order defined in COEFFICIENTS
            (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0, (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0,
            (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0, (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0,
            responseLossStaticSpecificCoefficient, (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0,
            (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0, (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0,
            (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0, (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0,
            (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0, (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0,
            (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0, (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0,
            (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0, (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0,
            (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0, (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0,
            responseInheritanceStaticSpecificCoefficient,
            (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0, (imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0,
        ],
        SA[ # order defined in COEFFICIENTS
            pathogenEstablishmentInteractionHostwideCoefficient, clearanceInteractionSpecificCoefficient, responseAcquisitionInteractionSpecificCoefficient,
            contactInteractionSpecificCoefficient, responseLossInteractionSpecificCoefficient,
            birthInteractionSpecificCoefficient, deathInteractionSpecificCoefficient, transitionInteractionSpecificCoefficient,
            receiveTransitionInteractionSpecificCoefficient, receiveContactInteractionSpecificCoefficient, intrahostFitnessInteractionSpecificCoefficient,
            generationsPerTransmissionInteractionSpecificCoefficient,
            mutationsPerGenerationInteractionSpecificCoefficient, recombinationsPerGenerationInteractionSpecificCoefficient, inoculumInteractionSpecificCoefficient,
            transmissionEfficiencyInteractionSpecificCoefficient, verticalTransmissionInteractionSpecificCoefficient, responseAcquisitionUponClearanceInteractionSpecificCoefficient,
            responseInheritanceInteractionSpecificCoefficient,
            (pat_g::String, imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0, (pat_g::String, imp_g::String, mat_g::String, hos_g::String, pop_id::String)->1.0,
        ],
        SA[ # order defined in COEFFICIENTS
            pathogenEstablishmentStaticHostwideCoefficient, clearanceStaticHostwideCoefficient, responseAcquisitionStaticHostwideCoefficient,
            contactStaticHostwideCoefficient, responseLossStaticHostwideCoefficient,
            birthStaticHostwideCoefficient, deathStaticHostwideCoefficient, transitionStaticHostwideCoefficient,
            receiveTransitionStaticHostwideCoefficient, receiveContactStaticHostwideCoefficient, intrahostFitnessStaticHostwideCoefficient,
            generationsPerTransmissionStaticHostwideCoefficient, mutationsPerGenerationStaticHostwideCoefficient, recombinationsPerGenerationStaticHostwideCoefficient, inoculumStaticHostwideCoefficient,
            transmissionEfficiencyStaticHostwideCoefficient, verticalTransmissionStaticHostwideCoefficient, responseAcquisitionUponClearanceStaticHostwideCoefficient,
            responseInheritanceStaticHostwideCoefficient, hostMutationsUponBirthStaticHostwideCoefficient, hostRecombinationsUponBirthStaticHostwideCoefficient,
        ],
        SA[ # order defined in COEFFICIENTS
            pathogenEstablishmentInteractionHostwideCoefficient, clearanceInteractionHostwideCoefficient, responseAcquisitionInteractionHostwideCoefficient,
            generationsPerTransmissionInteractionHostwideCoefficient, contactInteractionHostwideCoefficient, responseLossInteractionHostwideCoefficient,
            birthInteractionHostwideCoefficient, deathInteractionHostwideCoefficient, transitionInteractionHostwideCoefficient,
            receiveTransitionInteractionHostwideCoefficient, receiveContactInteractionHostwideCoefficient, intrahostFitnessInteractionHostwideCoefficient,
            mutationsPerGenerationInteractionHostwideCoefficient, recombinationsPerGenerationInteractionHostwideCoefficient, inoculumInteractionHostwideCoefficient,
            transmissionEfficiencyInteractionHostwideCoefficient, verticalTransmissionInteractionHostwideCoefficient, responseAcquisitionUponClearanceInteractionHostwideCoefficient,
            responseInheritanceInteractionHostwideCoefficient, hostMutationsUponBirthInteractionHostwideCoefficient, hostRecombinationsUponBirthInteractionHostwideCoefficient,
        ],
    )
end

function newHostType(
    id::String;
    template::HostType=DEFAULT_HOST_TYPE,
    num_loci::Union{Nothing,Int64}=nothing,
    possible_alleles::Union{Nothing,String}=nothing,

    # Each element takes seq argument, returns Float64
    pathogenEstablishmentCoefficient::Union{Nothing,Function}=nothing,
    clearanceCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionCoefficient::Union{Nothing,Function}=nothing,
    contactCoefficient::Union{Nothing,Function}=nothing,
    responseLossCoefficient::Union{Nothing,Function}=nothing,
    birthCoefficient::Union{Nothing,Function}=nothing,
    deathCoefficient::Union{Nothing,Function}=nothing,
    transitionCoefficient::Union{Nothing,Function}=nothing,
    receiveTransitionCoefficient::Union{Nothing,Function}=nothing,
    receiveContactCoefficient::Union{Nothing,Function}=nothing,
    intrahostFitnessCoefficient::Union{Nothing,Function}=nothing,
    generationsPerTransmissionCoefficient::Union{Nothing,Function}=nothing,
    mutationsPerGenerationCoefficient::Union{Nothing,Function}=nothing,
    recombinationsPerGenerationCoefficient::Union{Nothing,Function}=nothing,
    inoculumCoefficient::Union{Nothing,Function}=nothing,
    transmissionEfficiencyCoefficient::Union{Nothing,Function}=nothing,
    verticalTransmissionCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionUponClearanceCoefficient::Union{Nothing,Function}=nothing,
    responseInheritanceCoefficient::Union{Nothing,Function}=nothing,
    hostMutationsUponBirthCoefficient::Union{Nothing,Function}=nothing,
    hostRecombinationsUponBirthCoefficient::Union{Nothing,Function}=nothing,
)

    isnothing(num_loci) ? num_loci = template.num_loci : num_loci = num_loci
    isnothing(possible_alleles) ? possible_alleles = template.possible_alleles : possible_alleles = possible_alleles

    isnothing(pathogenEstablishmentCoefficient) ? pathogenEstablishmentCoefficient = template.coefficient_functions[PATHOGEN_ESTABLISHMENT] : pathogenEstablishmentCoefficient = pathogenEstablishmentCoefficient
    isnothing(clearanceCoefficient) ? clearanceCoefficient = template.coefficient_functions[CLEARANCE] : clearanceCoefficient = clearanceCoefficient
    isnothing(responseAcquisitionCoefficient) ? responseAcquisitionCoefficient = template.coefficient_functions[RESPONSE_ACQUISITION] : responseAcquisitionCoefficient = responseAcquisitionCoefficient
    isnothing(contactCoefficient) ? contactCoefficient = template.coefficient_functions[CONTACT] : contactCoefficient = contactCoefficient
    isnothing(responseLossCoefficient) ? responseLossCoefficient = template.coefficient_functions[RESPONSE_LOSS] : responseLossCoefficient = responseLossCoefficient
    isnothing(birthCoefficient) ? birthCoefficient = template.coefficient_functions[BIRTH] : birthCoefficient = birthCoefficient
    isnothing(deathCoefficient) ? deathCoefficient = template.coefficient_functions[DEATH] : deathCoefficient = deathCoefficient
    isnothing(transitionCoefficient) ? transitionCoefficient = template.coefficient_functions[TRANSITION] : transitionCoefficient = transitionCoefficient
    isnothing(receiveTransitionCoefficient) ? receiveTransitionCoefficient = template.coefficient_functions[RECEIVE_TRANSITION] : receiveTransitionCoefficient = receiveTransitionCoefficient
    isnothing(receiveContactCoefficient) ? receiveContactCoefficient = template.coefficient_functions[RECEIVE_CONTACT] : receiveContactCoefficient = receiveContactCoefficient
    isnothing(intrahostFitnessCoefficient) ? intrahostFitnessCoefficient = template.coefficient_functions[INTRAHOST_FITNESS] : intrahostFitnessCoefficient = intrahostFitnessCoefficient
    isnothing(generationsPerTransmissionCoefficient) ? generationsPerTransmissionCoefficient = template.coefficient_functions[GENERATIONS_PER_TRANSMISSION] : generationsPerTransmissionCoefficient = generationsPerTransmissionCoefficient
    isnothing(mutationsPerGenerationCoefficient) ? mutationsPerGenerationCoefficient = template.coefficient_functions[MUTATIONS_PER_GENERATION] : mutationsPerGenerationCoefficient = mutationsPerGenerationCoefficient
    isnothing(recombinationsPerGenerationCoefficient) ? recombinationsPerGenerationCoefficient = template.coefficient_functions[RECOMBINATIONS_PER_GENERATION] : recombinationsPerGenerationCoefficient = recombinationsPerGenerationCoefficient
    isnothing(inoculumCoefficient) ? inoculumCoefficient = template.coefficient_functions[INOCULUM] : inoculumCoefficient = inoculumCoefficient
    isnothing(transmissionEfficiencyCoefficient) ? transmissionEfficiencyCoefficient = template.coefficient_functions[TRANSMISSION_EFFICIENCY] : transmissionEfficiencyCoefficient = transmissionEfficiencyCoefficient
    isnothing(verticalTransmissionCoefficient) ? verticalTransmissionCoefficient = template.coefficient_functions[VERTICAL_TRANSMISSION] : verticalTransmissionCoefficient = verticalTransmissionCoefficient
    isnothing(responseAcquisitionUponClearanceCoefficient) ? responseAcquisitionUponClearanceCoefficient = template.coefficient_functions[RESPONSE_ACQUISITION_UPON_CLEARANCE] : responseAcquisitionUponClearanceCoefficient = responseAcquisitionUponClearanceCoefficient
    isnothing(responseInheritanceCoefficient) ? responseInheritanceCoefficient = template.coefficient_functions[RESPONSE_INHERITANCE] : responseInheritanceCoefficient = responseInheritanceCoefficient
    isnothing(hostMutationsUponBirthCoefficient) ? hostMutationsUponBirthCoefficient = template.coefficient_functions[HOST_MUTATIONS_UPON_BIRTH] : hostMutationsUponBirthCoefficient = hostMutationsUponBirthCoefficient
    isnothing(hostRecombinationsUponBirthCoefficient) ? hostRecombinationsUponBirthCoefficient = template.coefficient_functions[HOST_RECOMBINATIONS_UPON_BIRTH] : hostRecombinationsUponBirthCoefficient = hostRecombinationsUponBirthCoefficient

    return HostType(
        id,
        num_loci,
        possible_alleles,
        SA[ # order defined in COEFFICIENTS
            pathogenEstablishmentCoefficient, clearanceCoefficient, responseAcquisitionCoefficient,
            contactCoefficient, responseLossCoefficient,
            birthCoefficient, deathCoefficient, transitionCoefficient,
            receiveTransitionCoefficient, receiveContactCoefficient, intrahostFitnessCoefficient,
            generationsPerTransmissionCoefficient, mutationsPerGenerationCoefficient,
            recombinationsPerGenerationCoefficient, inoculumCoefficient,
            transmissionEfficiencyCoefficient, verticalTransmissionCoefficient, responseAcquisitionUponClearanceCoefficient,
            responseInheritanceCoefficient, hostMutationsUponBirthCoefficient, hostRecombinationsUponBirthCoefficient,
        ],
    )
end

function newPopulationType(
    id::String;
    template::PopulationType=DEFAULT_POPULATION_TYPE,
    constant_contact_density::Union{Nothing,Bool}=nothing,
    constant_transition_density::Union{Nothing,Bool}=nothing,
    host_sexual_reproduction::Union{Nothing,Bool}=nothing,
    hostSexualCompatibility::Union{Nothing,Function}=nothing,

    # Rate coefficients:
    pathogen_establishment_coefficient::Union{Nothing,Float64}=nothing,
    clearance_coefficient::Union{Nothing,Float64}=nothing,
    response_acquisition_coefficient::Union{Nothing,Float64}=nothing,
    contact_coefficient::Union{Nothing,Float64}=nothing,
    response_loss_coefficient::Union{Nothing,Float64}=nothing,
    birth_coefficient::Union{Nothing,Float64}=nothing,
    death_coefficient::Union{Nothing,Float64}=nothing,
    transition_coefficient::Union{Nothing,Float64}=nothing,
    intrahost_fitness_coefficient::Union{Nothing,Float64}=nothing,
    generations_per_transmission_coefficient::Union{Nothing,Float64}=nothing,
    mutations_per_generation_coefficient::Union{Nothing,Float64}=nothing,
    recombinations_per_generation_coefficient::Union{Nothing,Float64}=nothing,
    inoculum_coefficient::Union{Nothing,Float64}=nothing,
    transmission_efficiency_coefficient::Union{Nothing,Float64}=nothing,
    vertical_transmission_coefficient::Union{Nothing,Float64}=nothing,
    response_acquisition_upon_clearance_coefficient::Union{Nothing,Float64}=nothing,
    response_inheritance_coefficient::Union{Nothing,Float64}=nothing,
    host_mutations_upon_birth_coefficient::Union{Nothing,Float64}=nothing,
    host_recombinations_upon_birth_coefficient::Union{Nothing,Float64}=nothing,

    # Receive coefficients
    receive_transition_coefficient::Union{Nothing,Float64}=nothing,
    receive_contact_coefficient::Union{Nothing,Float64}=nothing,

    # Other functions
    pathogenFractions::Union{Nothing,Function}=nothing,
    # Takes Host entity, Population ID, and Population's weightedInteraction function,
    # returns vector with fractional representation of each pathogen present
    weightedInteractionPathogen::Union{Nothing,Function}=nothing,
    # Takes Pathogen entity, Host entity, Population ID, and event number;
    # returns aggregated response coefficient against that Pathogen for that event
    weightedInteractionResponse::Union{Nothing,Function}=nothing,
    # Takes Response entity, Host entity, Population ID, and event number;
    # returns aggregated response coefficient of that Response for that event
    weightedInteractionHostwide::Union{Nothing,Function}=nothing,
    # Takes Host entity, Population ID, and event number;
    # returns aggregated hostwide response coefficient based on all pairwise
    # interactions between Responses and Pathogens for that event

    developResponses::Union{Nothing,Function}=nothing,
    # takes in Pathogen, Host, Population ID,  population's Dict of Responses, population type's dictionary of ResponseTypes, and birth time as arguments, returns Response entities to be added
    # (this handles how many and which responses to choose when adding a response to a host)

    response_types::Union{Nothing,Dict{String,ResponseType}}=nothing,
)

    isnothing(constant_contact_density) ? constant_contact_density = template.constant_contact_density : constant_contact_density = constant_contact_density
    isnothing(constant_transition_density) ? constant_transition_density = template.constant_transition_density : constant_transition_density = constant_transition_density
    isnothing(host_sexual_reproduction) ? host_sexual_reproduction = template.host_sexual_reproduction : host_sexual_reproduction = host_sexual_reproduction

    isnothing(hostSexualCompatibility) ? hostSexualCompatibility = template.hostSexualCompatibility : hostSexualCompatibility = hostSexualCompatibility

    isnothing(pathogen_establishment_coefficient) ? pathogen_establishment_coefficient = template.base_coefficients[PATHOGEN_ESTABLISHMENT] : pathogen_establishment_coefficient = pathogen_establishment_coefficient
    isnothing(clearance_coefficient) ? clearance_coefficient = template.base_coefficients[CLEARANCE] : clearance_coefficient = clearance_coefficient
    isnothing(response_acquisition_coefficient) ? response_acquisition_coefficient = template.base_coefficients[RESPONSE_ACQUISITION] : response_acquisition_coefficient = response_acquisition_coefficient
    isnothing(contact_coefficient) ? contact_coefficient = template.base_coefficients[CONTACT] : contact_coefficient = contact_coefficient
    isnothing(response_loss_coefficient) ? response_loss_coefficient = template.base_coefficients[RESPONSE_LOSS] : response_loss_coefficient = response_loss_coefficient
    isnothing(birth_coefficient) ? birth_coefficient = template.base_coefficients[BIRTH] : birth_coefficient = birth_coefficient
    isnothing(death_coefficient) ? death_coefficient = template.base_coefficients[DEATH] : death_coefficient = death_coefficient
    isnothing(transition_coefficient) ? transition_coefficient = template.base_coefficients[TRANSITION] : transition_coefficient = transition_coefficient
    isnothing(receive_transition_coefficient) ? receive_transition_coefficient = template.base_coefficients[RECEIVE_TRANSITION] : receive_transition_coefficient = receive_transition_coefficient
    isnothing(receive_contact_coefficient) ? receive_contact_coefficient = template.base_coefficients[RECEIVE_CONTACT] : receive_contact_coefficient = receive_contact_coefficient
    isnothing(intrahost_fitness_coefficient) ? intrahost_fitness_coefficient = template.base_coefficients[INTRAHOST_FITNESS] : intrahost_fitness_coefficient = intrahost_fitness_coefficient
    isnothing(generations_per_transmission_coefficient) ? generations_per_transmission_coefficient = template.base_coefficients[GENERATIONS_PER_TRANSMISSION] : generations_per_transmission_coefficient = generations_per_transmission_coefficient
    isnothing(mutations_per_generation_coefficient) ? mutations_per_generation_coefficient = template.base_coefficients[MUTATIONS_PER_GENERATION] : mutations_per_generation_coefficient = mutations_per_generation_coefficient
    isnothing(recombinations_per_generation_coefficient) ? recombinations_per_generation_coefficient = template.base_coefficients[RECOMBINATIONS_PER_GENERATION] : recombinations_per_generation_coefficient = recombinations_per_generation_coefficient
    isnothing(inoculum_coefficient) ? inoculum_coefficient = template.base_coefficients[INOCULUM] : inoculum_coefficient = inoculum_coefficient
    isnothing(transmission_efficiency_coefficient) ? transmission_efficiency_coefficient = template.base_coefficients[TRANSMISSION_EFFICIENCY] : transmission_efficiency_coefficient = transmission_efficiency_coefficient
    isnothing(vertical_transmission_coefficient) ? vertical_transmission_coefficient = template.base_coefficients[VERTICAL_TRANSMISSION] : vertical_transmission_coefficient = vertical_transmission_coefficient
    isnothing(response_acquisition_upon_clearance_coefficient) ? response_acquisition_upon_clearance_coefficient = template.base_coefficients[RESPONSE_ACQUISITION_UPON_CLEARANCE] : response_acquisition_upon_clearance_coefficient = response_acquisition_upon_clearance_coefficient
    isnothing(response_inheritance_coefficient) ? response_inheritance_coefficient = template.base_coefficients[RESPONSE_INHERITANCE] : response_inheritance_coefficient = response_inheritance_coefficient
    isnothing(host_mutations_upon_birth_coefficient) ? host_mutations_upon_birth_coefficient = template.base_coefficients[HOST_MUTATIONS_UPON_BIRTH] : host_mutations_upon_birth_coefficient = host_mutations_upon_birth_coefficient
    isnothing(host_recombinations_upon_birth_coefficient) ? host_recombinations_upon_birth_coefficient = template.base_coefficients[HOST_RECOMBINATIONS_UPON_BIRTH] : host_recombinations_upon_birth_coefficient = host_recombinations_upon_birth_coefficient

    isnothing(pathogenFractions) ? pathogenFractions = template.pathogenFractions : pathogenFractions = pathogenFractions
    isnothing(weightedInteractionPathogen) ? weightedInteractionPathogen = template.weightedInteractionPathogen : weightedInteractionPathogen = weightedInteractionPathogen
    isnothing(weightedInteractionResponse) ? weightedInteractionResponse = template.weightedInteractionResponse : weightedInteractionResponse = weightedInteractionResponse
    isnothing(weightedInteractionHostwide) ? weightedInteractionHostwide = template.weightedInteractionHostwide : weightedInteractionHostwide = weightedInteractionHostwide
    isnothing(developResponses) ? developResponses = template.developResponses : developResponses = developResponses

    isnothing(response_types) ? response_types = template.response_types : response_types = response_types

    return PopulationType(
        id,
        constant_contact_density,
        constant_transition_density,
        host_sexual_reproduction,
        hostSexualCompatibility,
        SA[ # order defined in COEFFICIENTS
            pathogen_establishment_coefficient, clearance_coefficient, response_acquisition_coefficient,
            contact_coefficient, response_loss_coefficient,
            birth_coefficient, death_coefficient, transition_coefficient,
            receive_transition_coefficient, receive_contact_coefficient, intrahost_fitness_coefficient,
            generations_per_transmission_coefficient, mutations_per_generation_coefficient,
            recombinations_per_generation_coefficient, inoculum_coefficient,
            transmission_efficiency_coefficient, vertical_transmission_coefficient, response_acquisition_upon_clearance_coefficient,
            response_inheritance_coefficient, host_mutations_upon_birth_coefficient, host_recombinations_upon_birth_coefficient,
        ],
        pathogenFractions,
        weightedInteractionPathogen,
        weightedInteractionResponse,
        weightedInteractionHostwide,
        developResponses,
        response_types,
    )
end
