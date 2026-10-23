# jOpqua Changelog

## Known Issues:

### Known bugs:
- Inter-population contact rates are higher than intra-population contact rates
- Strange increase in infected population towards end of simulation in 
`pathogen_evolution_recombination_test.jl`
- Host birth recombinations are broken and use deprecated parameters

### Known performance issues:
- Garbage collection time is higher than expected for 100k host simulations, suspected type
instabilities in weights/immunity changes
- Births/deaths result in significant slowdown due to dynamically resizing (thus re-declaring)
host weight matrices, as well as garbage collection (probably associated to the former)

### Mechanistic/biological issues and limitations:
- [New paper](https://www.biorxiv.org/content/10.1101/2023.11.19.567585v3) shows transmission 
bottlenecking is more complex than thought: an initial infection period captures a high
diversity of genetic variability from the donor, and a subsequent intrahost bottleneck
presumably due to overcoming first line host defenses results in a single lineage founder
event for the infection. Intrahost selection can act at that bottleneck. We could build this
into the transmission event function by sampling based on fitness within the new host.
- Within-host dynamics are currently approximated as stochastic Poisson processes with variable 
rates or through the use of compartment models (e.g. Exposed and Infected compartments in an
SEIR-type model) by using separate `Population` entities for each compartment and setting 
transition rates between them. However, detailed `Pathogen` and `Response` population changes 
within individual hosts are not currently captured.
- Similarly, within-host lineage establishment events in which multiple mutations and/or 
recombinations co-occur do not take into account the fitness of intermediate genotypes and 
may be inaccurate.

### TODO:
- Use `population_id` `String` as function argument everywhere needed where 
we can't use `Population` as an argument
- Reexamine clearance rate time correction
- Change transmission event function to add receiver-side bottleneck option
- Change ancestor search so that you can choose to get only a fraction of ancestors 
(improves the function's runtime for simulations with long timelines)
- Try explicitly-declared max host population size in instance of population change,
dynamically resize if max population exceeded

Debug the following:
- Mutant establishment
- Recombinant establishment
- Diploid `Host` genomes/homologous chromosome separators
- `Host` birth events: vertical transmission, `Response` inheritance, mutation 
and recombination upon birth
- Interventions

## Changelog:

## 22 Oct 2026
- Fixed floating point error in total rate sum for populations causes error in 
some simulation runs, e.g. `metapopulations.jl`
- Added `recalculation_counters` to model to avoid floating point error
- Changed `randChoose()` zero-sum error handling for function typing, to avoid 
returning `NaN` and messing up `Int64` type
- Fixed simple bug in `propagateWeightsOnRemoveHost()` where there was an extra
`model.`

## 20 Oct 2026
Debugging floating point! WIP

## 4 Jun 2026
- Changed argument order for immunity sequence functions
- Played around with immunity examples
- Made `Pathogen` sequence coefficient functions take `Population` ID as 
additional argument
- Made `Response` sequence coefficient functions take `Population` ID as 
additional argument
- Made `Host` sequence coefficient functions take `Population` ID as 
additional argument

Considered making `Pathogen` sequence coefficient functions take `Host` genome 
as an additional argument. However, decided against it to preserve 
`Pathogen`-level coefficients as intrinsic. If we want effects specific to 
different `Host` genomes, we use `Response` entities.

Considered moving `pathogenFractions`, `weightedInteractionPathogen`, 
`weightedInteractionResponse`, `weightedInteractionHostwide`, 
`developResponses`, and `response_types` into `HostType` fields with 
additional `Population` ID as argument, but not technically possible because 
of circular definition in `HostType` and `Host`

Added `test.jl` to test out possible algorithms to handle clearance events; 
this is just a scratchpad for experimenting.

## 28 May 2026
- Fixed argument bug in `newResponseType()`

## 27 May 2026
- Created `choosePathogenFromPopulationFraction()` as a function in `choice.jl` to 
remove `randChoose()` function calls from `events.jl`
- Renamed `LINEAGE_ESTABLISHMENT` to `PATHOGEN_ESTABLISHMENT` and all associated 
parameters

Another overhaul of clearances:

- Changed `hostWeights()` to compute `Host`-level `CLEARANCE` rate to be sum of 
`Pathogen`-level `CLEARANCE` rates and thus avoid long coinfection duration
- Added `Host` property `type_clearance_numbers` to store a dictionary 
where each `PathogenType` is mapped to the sum of total clearance weights 
of all `Pathogens` with and without the extinction rate due to drift and 
the count of `Pathogens` of each `PathogenType`
- Changed computation of `CLEARANCE` within `pathogenWeights()` to store the 
above information for each `PathogenType` and also add the rate of extinction 
due to drift
- Changed `hostWeightsPathogen()` to compute `CLEARANCE` weights the same as 
other weights
- Changed `clearPathogen()` to check whether a whole `PathogenType` is 
to be cleared and execute it if so
- Added `pathogenClearanceThroughDriftWeight()` to compute the clearance weight 
due to drift for a `Pathogen`, currently just a placeholder that returns 0.0

This breaks the assumption of "splitting the immune response" from earlier; 
cross-immunity should handle this. Considered whether to implement cross-immunity 
check at clearance, but decided the best way to handle this is to sum clearance 
rates of different pathogens within a `Host` (representative of their smaller 
within-host population fraction and increased extinction likelihood through 
drift/competitive exclusion). What should definitely not sum and remain proportional 
to pathogen population fraction is the `RESPONSE_ACQUISITION` rate. When a 
`Response` is acquired, then we have cross-immunity and the clearing rates of 
different pathogens within a `Host` may become correlated.

I'm not convinced by this, though. It's worth pointing out that there are different 
types of clearance: 
- Clearance of all `Pathogens` of the same `PathogenType` within a `Host`, which may 
be due to non-specific, innate immunity or general changes in `Host` environment
and would be computed by weighting `Pathogen` `CLEARANCE` coefficients by the 
within-host population fraction (weighted average of infecting `Pathogens`) of the 
same type (which itself would involve one fraction count per `PathogenType` within
each `Host`; pretty impractical)
- Extinction of a specific `Pathogen` due to drift and exclusion by other `Pathogens` 
of the same `PathogenType` within the `Host`, which is a function of the total
niche size for that `PathogenType` and may be affected by the same non-specific, 
innate immunity and general changes in `Host` environment mentioned above, and 
would have a total per `Host` rate equal to the sum of rates for each `Pathogen` 
within it
- Clearance of a specific `Pathogen` due to acquired immunity effects specific to 
that `Pathogen`, determined by `Response` entities within the `Host`, which also
would have a total per `Host` rate equal to the sum of rates for each `Pathogen` 
within it

We can choose to handle them all as separate events, but that's innacurate because 
they actually act in concert with each other: specific immune effects can affect 
the population fraction of a `Pathogen` lineage and increase the likelihood that 
it's lost to drift, for example. This could be accounted for through hostwide 
coefficients in `Response` entities, though. It would imply adding two new types 
of coefficients and events to the simulation, though. Additionally, the multiple 
within-host fraction counts sounds like a terrible idea, especially when 
considering different degrees of niche-sharing between different `PathogenTypes`.

Instead, since all this was prompted specifically because of the problem of 
complex infections having lower clearance rates than single pathogen infections, 
I'm choosing to use the following approximation:

Clearances are considered to be a single type of event that selects a specific 
`Pathogen` in a `Host`. The total per `Host` clearance weight is equal to the 
straight sum of `CLEARANCE` weights for each `Pathogen` within it (not 
weighted by population fraction because fraction could affect the weight in 
different ways, either increasing `CLEARANCE` weight due to low populations or 
decreasing it due to low immune profile). Each `Pathogen`'s clearance weight 
is equal to the sum of two components: 
1. The product of its base clearance coefficient with all `Pathogen` and 
`Response` effect factors present within the `Host`
2. The weight (or rate) of extinction due to drift, computed from population 
genetics based on population size, generation time, and fitness (equal to 
zero in cases of single infections)

Each `Host` computes a number reflecting the mean of all `Pathogen` clearance 
weights not associated with drift (and therefore not just an consequence of 
the basic notion of coinfection, i.e. the mean of all weights from 
Point 1 above) for each `PathogenType`. This is what the total clearance 
weight (that is, leading to the loss of all pathogens of that `PathogenType`)
would be if we treated them as a single population of pathogens, which makes 
sense as an assumption since they're all presumably sharing similar resources 
and affected by the same non-specific immune responses at least to some degree.
In other words, to solve the discrepancy in the total clearance times of 
simple and complex infections, this is what the total clearance weight should 
be. That number is then divided by the sum of all `Pathogen` clearance 
weights for the given `PathogenType` (including the part of each weight 
derived from drift). This quotient is then stored for each `PathogenType` 
present in the `Host`, or even just for each `Pathogen`, even if it's 
repetitive information. Call it `type_clearance_probability[pathogen]`.

To make the timelines of complex and simple infections align, each `CLEARANCE` 
event begins with a check of probability `type_clearance_probability[pathogen]`
(with `type` being the chosen `Pathogen`'s `PathogenType`). If the 
check succeeds, all `Pathogens` of the chosen `PathogenType` within the `Host` 
are cleared. Otherwise, only the specific `Pathogen` chosen in the event is 
cleared. This guarantees that at every point in time, the net, total clearance 
rate for any `PathogenType` within the `Host` is equal to the mean clearance 
rate for that `PathogenType`. Despite this, `Pathogen`-specific variation 
due to intrinsic clearance coefficient, `Pathogen` and `Response` effects, and 
extinction through drift is preserved.

Admittedly, this approach may result in unexpectedly fast clearance times when 
two `Pathogens` of the same `PathogenType` but vastly different `CLEARANCE` 
weights infect the same `Host` simultaneously. One could expect the total 
time to be more strongly determined by the longer-lived `Pathogen`. However,
the assumption here is that the `Pathogen` with the higher `CLEARANCE` weight 
is consuming more of the `Host`'s resources shared with the longer-lived 
`Pathogen` and triggering a stronger innate response that also affects the 
longer-lived `Pathogen`. A helpful analogy is to consider `DEATH` events, 
which do not have a wait time determined by the longest-lived `Pathogen`
but by the fraction-weighted average of all infecting `Pathogens` (in the 
case of `CLEARANCE` events, no fraction weighting is applied as explained
above, however).

I'm also considering alternative approaches though. This will update soon.

## 25 May 2026
- Fixed bugs in `recombinantSequences()` that affected proper generation of recombinant 
sequences
- Fixed syntax errors in `establishRecombinant!()` (later removed)
- Fixed Julia syntax error for `String` concatenation in `recombinantSequences()` and 
fixed bug with `String`-based `MVector`s not being able to be index-set

Made lineage establishment to be an event in which recombination and a certain number of 
mutations may happen. This required a medium overhaul:
- Changed `MUTATIONS_UPON_INFECTION` to `MUTATIONS_PER_GENERATION` and 
`RECOMBINATIONS_UPON_INFECTION` to `RECOMBINATIONS_PER_GENERATION`
- Removed `RECOMBINANT_ESTABLISHMENT` as a coefficient type
- Added `GENERATIONS_PER_TRANSMISSION` as a `Pathogen`-specific non-sampling variable 
within the model coefficients
- Changed `MUTANT_ESTABLISHMENT` to `LINEAGE_ESTABLISHMENT`
- Changed `DEFAULT_POPULATION_TYPE` to have a starting base coefficient of `1.0` for 
`GENERATIONS_PER_TRANSMISSION`
- Changed `pathogenWeights!()` to adjust the `LINEAGE_ESTABLISHMENT` weight in cases in 
which there is a single pathogen within the host and effective recombination is not 
possible
- Changed `mutantPathogen!()` and `recombinantPathogens!()` to accept parameters with 
the expected number of mutations and recombinations per generation
- Created `establishLineage!()` based on `establishMutant!()` to handle lineage 
establishment events from either mutation or recombination
- Removed `establishRecombinant!()`
- Changed `hostContact!()` to correctly use syntax of new `mutantPathogen!()` and 
`recombinantPathogens!()` functions, and to account for the number of generations 
elapsed before transmission

Originally, I considered calculating the number of mutations and recombinations based on 
the number of generations elapsed before establishment, which itself can be obtained 
based on relative rates of lineage establishment and infection (or rather the inverse, 
time to establishment and time to infection) along with the number of generations per 
transmission. However, we are choosing a different approach, at least for lineage 
establishment: if we establish a new lineage, we assume the new genotype is the center 
of a mutant cloud, so we only consider lineage establishment with the mutations and 
recombinations accumulated in a single generation at a time.

We also assume that the rate of lineage establishment comprehends both the rate of 
establishment of a new lineage and the rate of lineage emergence through both mutation 
and recombination in a host with a perfectly diverse pathogen population (i.e., one 
within which any recombination is effective recombination). This allows us to subtract
the fraction of lineage establishment that occurs from recombination in hosts with 
only a single pathogen.

### 24 May 2026
Major update solving multiple important bugs related to transitions and contacts. 

Inter-population contacts and transitions considered debugged!

Changed calculation of `RECEIVE_TRANSITION` weight in `population_weights_receive` to use 
mean of host weights instead of sum, given that receiving transitions doesn't scale with 
the number of hosts. This requires changes:
- Change in `newPopulation!()` to set the starting receive transition weight of a new 
empty population to base coefficient value
- Changed `calculateWeightReceive!()` to set the receive transition weight of an empty
population to the base coefficient value
- Changed `calculateWeight!()` to include `population_contact_weights_receive_sums` in 
`population_weights` calculation for `TRANSITION`
- Changed `Host`-level `propagateWeightChanges!()` to update `host_weights` with 
the `change` value before `change` is itself changed by more factors in cases of 
contact or transition events to modify `host_weights` without including 
`population_transition_weights_receive_sums` (also 
`population_contact_weights_receive_sums`)
- Changed default starting coefficient for receiving transitions to be `1.0` in 
`DEFAULT_POPULATION_TYPE`
- Changed default value of `constant_transition_density` to `true` in 
`DEFAULT_POPULATION_TYPE`; I'm not sure why it was `false` to begin with
- Changed `updatePopulationTransitionWeightReceiveMatrix!()` to not update 
`population_transition_weights_receive_sums`, since that is handled in the internal
`calculatePopulationTransitionWeightReceiveMatrix!()` call
- Changed `updatePopulationContactWeightReceiveMatrix!()` to not update 
`population_contact_weights_receive_sums`, since that is handled in the internal
`calculatePopulationContactWeightReceiveMatrix!()` call
- Overhauled `propagateWeightsOnAddHost!()` to check for constant density in `if` 
statements when deciding to adjust for transition and contact event denominators,
propagating weights with the correct factors in those cases, temporarily adjusting 
the transition rates to zero if an empty population is getting a new host, and using 
`host_weights` instead of directly calculating and triggering weight propagations
from within the function (this was incorrect because it wasn't accounting for host
coefficients modified by host genome, pathogens, or responses)
- Changed `metapopulation.jl` to illustrate inter-population contact and transition 
behavior

Remember that `population_weights_receive` doesn't divide by population size using the 
`constant_transition_density` parameter, that is only included in the 
`population_contact_weights_receive` matrix.

Also:
- Added `@views` to `host_weights` and `host_weights_receive` slice accessors in 
`removeHostFromPopulation!()` and `calculateWeight!()`
- Removed `START_COEFFICIENT` and the two references to it (they were messing up rate 
calculations)
- Changed `addHostToPopulation!()` to use `propagateWeightsOnAddHost!()` to update 
population weights
- Changed `addHostsToPopulation!()` to calculate new `Host` weights once and use them 
for all added hosts
- Created `propagateWeightsOnRemoveHost!()` to update weights when a host is removed
- Changed `removeHostFromPopulation!()` to use `propagateWeightsOnRemoveHost!()` to 
update population weights

Later:
- Added back recalculation counters for receive weights
- Removed redundant resetting of recalculation counters in `propagateWeightChanges!()` and 
`propagateWeightReceiveChanges!()`
- Changed parameters of `pathogen_evolution.jl` example to get nice selective sweeps again,
this time there's some interesting clonal interference too

### 16 May 2026
- Modified parameters of `pathogen_evolution.jl` example to get nice selective sweeps again

Multiple bug fixes for multi-population simulations:
- Fixed bug in `newPopulation!()` when expanding population contact and transition matrix 
while adding new populations
- Fixed bug in `saveCompartments()` when constructing compartment matrix without assigning
times to all populations
- Fixed bug in `plotCompartments()` where time vector was wrong size (deleted repeats from
multiple populations)
- Fixed bug in `propagateWeightsOnAddHost!()` which miscalculated the change in 
`population_contact_weights_receive_sums` when multiple populations are present

Transition and contact:
- Changed syntax of `setPopulationTransitionCoefficient!()` and 
`setPopulationContactCoefficient!()` to use population objects rather than population IDs,
aiming to make the API consistently use objects directly rather than IDs
- Changed default starting coefficient for receiving transitions to be `1.0` in 
`START_COEFFICIENTS`
- Fixed bug in `simulate!()` to be able to deal with empty populations 

### 14 May 2026
- Modified clearance weight calculation to be the maximum clearance weight of all pathogens
within the host, accounting for responses, rather than the sum. This solves the issue where
coinfections of antigenically similar strains have overestimated clearance rates. 
- Choosing not to implement a system that checks for cross-immunity when acquiring immunity 
upon clearance; this would imply lots of random number generation to evaluate each present 
strain as well as some mapping of cross-immunity to probability. We can just let the 
simulation figure it out based on the recomputed clearance rates.

Together, this results in slightly lower net clearance rates of individual strains in 
coinfections with other distantly related strains, but that is a reasonable assumption 
representing "splitting the immune system's total resources".

Original issue description:
- Conceptual issue: Current way of handling intrahost pathogen populations as fractions causes
clearance math to be off for complex infections (coinfections). Imagine two strains with equal
fitness and 100% cross-immunity to each other coinfecting an individual. For starters, in the
current model, coinfection clearance is now dependent on two independent Poisson processes,
an incorrect assumption to begin with since the second clearance event is much more likely
once the first occurs. This results in the clearance time of the coinfection in the model
now following a Gamma distribution rather than an exponential one, which should not be the
case, but might be difficult to work around. However, the issue is not only in the
distribution of clearance times, but also in their mean. In the current model, the coinfected
individual's clearance rate is twice that of that of the individual with the single infection.
This is correct as far as the mean rate goes, since each strain has half the population it
would have if it were infecting an individual with no other strains present. However, once the
first strain is cleared, the clearance rate of the second strain stays the same, so the total
clearance rate for that individual host is now half of what it was before the first clearance
event. This results in a net increase in mean clearance rate for the coinfection with respect
to the single infection--a clear mistake. To fix this, we can do the following two things:
(1) Modify clearance weight calculation to be the maximum clearance weight of all pathogens
within the host, accounting for responses, rather than the sum (this results in slightly lower
net clearance rates of individual strains in coinfections with other distantly related
strains, but that is a reasonable assumption representing "splitting the immune system's total
resources")
(2) Modify the clearance event to not only check whether it induces acquisition of a response,
but if it does, to also to result in clearance of other coinfecting pathogens if they cannot
escape the new immune responses (a probability determined by cross-immunity to those acquired
responses)... Maybe unnecessary, let subsequent clearances happen based on new clearance rates

### 10 May 2026
Back at it.
- Removed redundant parameters: `PathogenType` coefficient functions that are 
specific to `Response` events/nonsampling coefficients or to `Host` 
events/nonsampling coefficients, as well as `ResponseType` coefficient functions 
that are specific to `Pathogen`. Left function array structure intact by having the 
setup function for `PathogenType` and `ResponseType` fill in the redundant 
coefficient functions with functions that map to `1.0`.
- Changed parameters of `pathogen_evolution` simulation to show sweeps and clonal 
interference

### 21 Jul 2025
- Remove external setting of Flexle sums (reduntant with Flexle `1.0.3`)

### 15 Jul 2025
- In-house Hamming distance function for both (1) performance and (2) generalizability 
(e.g. nucleotide/amino acid substitution matrix in place of simple pairwise inequality; 
not currently implemented) (CLM)

### 14 Jul 2025
- Upgrade Flexle to `1.0.3`, fixes an uncommon bug when updating some weights

### 16 Jun 2025
Float precision errors in rate calculation accumulate over time. To solve this, I
added a full rate recalculation procedure that runs every certain number of events
to reset error (including for the `.sum` field on Flexle samplers, the `contact_sum` and
`transition_sum` fields in `Population` structs, and the sums in the `Model` struct)

- Changed weight propagation functions to check if change is zero and only do stuff if it
isn't
- Renamed Flexle sampler vector fields as plural samplers
- Added a vector field in `Population` to track the number of events since the last
weight sum recalculation
- Changed weight propagation to recalculate weight sums at both `Population` and `Model`
levels whenever a threshold is exceeded in the recalculation counter
- Added said threshold as a constant, `RECALCULATION_EVENT_THRESHOLD`

Beauutiful selective sweep succession in `pathogen_evolution` example.

### 15 Jun 2025
Finally tracked down a bug involving floating point errors in total sum calculation causing
clearances to be attempted when there was nothing to clear, particularly in coinfections.
To do this, I:
- Changed the behavior of `becomesZero` to check for fractional change in value rather
than absolute difference
- Reset Flexle sampler sums to freshly calculated total whenever a recalculation is
triggered by `becomesZero`

This then begs the question on whether Flexle should internally be checking for large changes
using `becomesZero` internally and recalculating on its own.

### 12 Jun 2025
- Added approximately zero check in `updatePopulationTransitionWeightReceiveMatrix!`
- Added coefficient recalculation functions to `weights.jl`, triggered when a coefficient sum
equals approximately zero (and didn't before)
- Added `becomesZero` as a `utils.jl` function for the above purpose

### 9 Jun 2025
Births and deaths now work.
- Bug fix in death handling within `removeHostFromPopulation!`
- Change `simulate` to not crash when host population is zero
- Change `saveHistory` to not crash when host population is zero
- Syntax errors and multiple other bug fixes and updates in `birth!`
- Removed `UMAP` requirement

### 7 Jun 2025
- Small change in `developResponse` structure, returning a single `Response` instead of a
vector

### 6 Jun 2025
- Small change in simulation algorithm to avoid running a final event if event time is over
max simulation time

### 5 Jun 2025
- Added `ancestors` function to data analysis package to pull out all ancestors from a list
of `Pathogen` genomes in a `Population`
- Multiple edits to flu example, code for genetic and antigenic maps might be cleaned up and
standardized into a `plots` function

I believe Hill function does now behave as expected based on all-none immunity (as of
previous commit).

### 24 May 2025
- Fixed bug in `weightedInteractionResponse` parameter types
- Fixed bug in `weightedInteractionHostwideProduct`, called the wrong function (specific
instead of hostwide)

### 15 May 2025
- Changed parameter in `pathogen_evolution.jl` example to the new correct hostwide parameter
- Removed type instabilities in `immunity.jl` functions

### 14 May 2025
- Changed order of loading files to avoid import conflicts
This is actually crucial for performance! Without it, the new immunity changes led to a ~3.8X
slowdown for the `pathogen_evolution.jl` with 10000 hosts

### 13 May 2025
- Continued changing nomenclature of specific/hostwide and static/interaction coefficients
- Changed weight computation to correctly account for all kinds of coefficients
- Added new parameters for functions controlling `Pathogen`-`Response` interaction weighting
for a specific `Pathogen`, a specific `Response`, and a hostwide function
(`weightedInteractionPathogen`, `weightedInteractionResponse`, and `weightedInteractionHostwide`,
respectively)

### 11 May 2025
- Corrected bug in `weightedInteractionWinnerTakesAll`
- Renamed "specific" (interaction-specific) coefficients as "interaction" coefficients
- Renamed "normal" coefficients are "specific" coefficients (specific to a pathogen or
response within a host, as opposed to "hostwide" coefficients)

There are now clearly seven types of coefficients modifying `Population`-level parameters:
1. `Pathogen` specific coefficients
2. `Pathogen` hostwide coefficients
3. `Response` specific static coefficients
4. `Response` specific specific coefficients
5. `Response` hostwide static coefficients
6. `Response` hostwide specific coefficients
7. `Host` (hostwide) coefficients
All parameters and coefficients are specified in the respective "type" structs:
`PathogenType`, `ResponseType`, `HostType`, `PopulationType`.

### 8 May 2025
- Added legend omission to plots
- Added Hill function to utils

### 25 April 2025
- Add `Flexle.sample` alias to resolve `sample` function name conflict with `StatsBase`

### 16 April 2025
- Removed local Flexle files
- Moved Flexle examples to outer directory
- Updated Flexle requirement to `1.0.2`
- Added code to use imported Flexle package (including `Flexle.sample` to resolve function
naming conflicts with `sample`)

### 12 April 2025
- Incorporated Flexle fixes implemented by CLM, results in 63.6% improvement in runtime
over `randChoose` linear search for `pathogen_evolution.jl` with 10,000 hosts; was even
able to run that simulation with 10^6 (!!!) hosts in 1034 seconds!

### 11 April 2025
- Moved `flexleSamplers` methods from `Flexle` into `utils.jl`

### 10 April 2025
- `Flexle` sampling has been incorporated into `Host` sampling within a `Population`;
`Flexle` results in performance decrease compared to `randChoose` due to possible
inefficiencies in the calling of `maxLevelWeight` (PCR)

### 29 March 2025
- Remove `Flexle` source code, now in [separate package](https://github.com/connormurphy798/Flexle.jl)
included in `jOpqua` dependencies (CLM)

### 27 March 2025
- Minor bug fix (`Flexle.index_positions` not updating when 0-weight element pushed) (CLM)

### 26 March 2025
- Precompute log2 of `Flexle` max level upper bound for faster indexing (CLM)

### 25 March 2025
- `Flexle.setindex!` avoid adding to/removing from level when new and old weights belong
in same level (CLM)

### 22 March 2025
- Reorganize `Flexle` into its own module (CLM)
- Refactor API to be more idiomatically "Julia", i.e. rename functions where appropriate
to be methods of existing functions like `push!` or `deleteat!` (CLM)
- Write `Flexle` API doc (CLM)

### 21 March 2025
- New strategy for recording `index_positions` as `Vector` at `FlexleSampler` level
(as opposed to `Dict` at `FlexLevel` level), ~30-50x speedup to addition/removal of
new weights (CLM)
- Revert to old `rejectionSample` technique; yet-unidentified bug in previous
strategy caused sample distribution to be wrong

### 19 March 2025
- Added (+clarified existing) docstrings for most `flexle.jl` functions (CLM)

### 7 March 2025
- Incorporated static and interaction hostwide coefficients of `Pathogen` and
`Response` entities into weight and coefficient propagation from the
`Pathogen` and `Response` level to the `Host` level
- Removed redundant `transmissionEfficiency` functions
- Harnessed calculated coefficients ratehr than functions in `pathogenWeights!`
- Fixed bugs regarding multiple mutations in `mutantSequence!`

### 6 March 2025
- Rename `weightedResponse` to `weightedInteraction`
- Fixed some bugs in `birth!`, changed host contact function to determine inoculum
from `INOCULUM` coeefficient by default but from `VERTICAL_TRANSMISSION`
coefficient during vertical transmission instead
- Make `Host` sequence coefficients affect sampling event coefficient propagation
- Futher improvements to `removeFromFlexleSampler!` dict update via streamlined
algorithm and fast hashing keys (CLM)

### 5 March 2025
- Add response acquisition upon clearance as preferred alternative to response
acquisition during infection (as in mutations upon infection vs. mutation
establishment)
- Transfer vertical transmission, response inheritance, host mutation, and host
recombination to nonsampling variable/coefficient/function vectors
- Improvement to `removeFromFlexleSampler!` runtime via non-brute force algorithm
(~33% speedup) (CLM)
- Remove unnecessary comments, correct long description in README from last
weekend

### 4 March 2025
- Renamed `infectionProbability` to `transmissionEfficiency`
- Transfer transmission efficiency and infection probability to nonsampling
variable/coefficient/function vectors
- Fix `Flexle` bug where `index_positions` was not updated on element removal; current
fix disastrously slow, needs speedup (CLM)

### 3 March 2025
- Added nonsampling coefficient vector to `Host` and some machinery to update it
- Create `HostType` struct to house `Host`-specific parameters
- Changed `birth!` to carry out a contact event if vertical transmission happens
- Transfer mutations upon infection, recombinations upon infection, and inoculum
to nonsampling variable/coefficient/function vectors
- Additional x2 weight update performance improvement via: `Float64` bit manipulation
to calculate log bounds + floor log2; and iterative approach to `levelIndex` (CLM)

### 2 March 2025
- Corrected a couple references from to `NUM_COEFFICIENTS` for
`INTRAHOST_FITNESS`
- Added new coefficients to hierarchy matrices in preparation of overhaul to
non-sampling variables and sub-events

### 1 March 2025
Some conceptual simulation structure advances:

An integral part of jOpqua simulations are the myriad of variables that affect
the course of the simulation based on the sequences of the `Pathogen`, `Host`
and `Response` events involved. These variables include numerical values such
as coefficients, weights, rates, probabilities, and probability distribution
parameters, as well as the sequence-evaluating functions that define those
parameters. These values and their respective functions can be classified into
two categories:

1. Those that affect event rates and determine sampling weights of entities
participating in those events ("sampling variables")
2. Those that determine the probabilities of "conditional sub-events" that may
be triggered when a main simulation event happens ("non-sampling variables")

The first category of values and functions must be fully computed at all times
in order to accurately determine the next event to be carried out and the
entities participating in it. This is why all entities in the hierarchy contain
matrices in which `Pathogen`, `Response`, `Host`, and `Population` weights are
ordered in a way such that they can be used to sample the corresponding entity.
Changes in these values are propagated up the hierarchy according to propagation
functions that incorporate relevant factors to be multiplied as you move up the
hierarchy.

The second category of variables do not modify the rates of events or the weights
used to sample entities in the simulation. Therefore, there is no need to
calculate the weights of individual entities for sampling. Furthermore, the exact
values of the probabilities and distribution parameters involved do not require
being fully computed and can be calculated individually on the fly as they become
needed. This may be more efficient than storing all the computed coefficients and
probabilities/distribution parameters, and then propagating changes to all
affected coefficients and probabilities whenever an entity changes, as done with
sampling variables. However, this is only true as long as the rate at which events
involving the same entities are sampled multiple times with no change in the
entities is very low (a reasonable assumption most of the time, I would argue).
Additionally, it requires less memory and is more straightforward to implement.
Thus, these coefficient functions will be stored in static vectors (could be
dictionaries since the order doesn't matter, but static vectors should
theoretically improve performance) and a general, full calculation function
computes the resulting probability or distribution parameter based on all the
involved coefficients and values. For each non-sampling probability or
distribution parameter computed, these coefficients and values include:

- Base probability or distribution parameter value set at the `PopulationType`
level (just one factor stored in `PopulationType`)
- Static, specific coefficients calculated based on a single sequence at the
`Pathogen`, `Response`, and `Host` level using functions defined in
`PathogenType`, `ResponseType`, and `HostType`, respectively, and involving
only the specific entities (including `Host` container) involved in the event
(one factor per entity involved, stored in the corresponding `Pathogen`,
`Response`, or `Host`)
- Static, hostwide coefficients calculated based on a single sequence at the
`Pathogen`, `Response`, and `Host` level using functions defined in
`PathogenType`, `ResponseType`, and `HostType`, respectively (one factor,
stored at the `Host` or `Population` levels for `Pathogen`/`Response` and
`Host` level aggregation, respectively, aggregated from the list of
coefficients of each entity present in the same `Host` based on
`pathogenFraction` for pathogens or `responseAggregate`--which by default is
the total product--for responses, or `Population` based on
`hostAggregate`--the default of which is just `1.0`--, respectively)
- Interaction-based, specific coefficient calculated based on four sequences
(one for `Pathogen`, two for `Response`, and one for `Host`) using a function
defined in `ResponseType` (one factor, computed on the fly, aggregated from
the list of coefficients of each `Pathogen`-`Response` combination including
either the specific `Pathogen` or `Response` involved in the non-sampling
sub-event present within the `Host`, aggregated within a `Host` based on
`weightedInteractionSpecific`)
- Interaction-based, hostwide coefficient calculated based on four sequences
(one for `Pathogen`, two for `Response`, and one for `Host`) using a function
defined in `ResponseType` (one factor, stored within the `Host`, aggregated
from the list of coefficients of each `Pathogen`-`Response` combination present
within the `Host`, aggregated within a `Host` based on
`weightedInteractionGeneral`)

The list of non-sampling coefficient functions should include (as far as I have
thought), in alphabetical order:

- `HOST_MUTATIONS_UPON_BIRTH` (modifies the mean of a Poisson random variable)
- `HOST_RECOMBINATIONS_UPON_BIRTH` (modifies the mean of a Poisson random variable)
- `INOCULUM` (modifies the mean of a Poisson random variable)
- `MUTATIONS_UPON_INFECTION` (modifies the mean of a Poisson random
variable)
- `RECOMBINATIONS_UPON_INFECTION` (modifies the mean of a Poisson
random variable)
- `RESPONSE_ACQUISITION_UPON_CLEARANCE` (modifies a probability)
- `RESPONSE_INHERITANCE` (modifies a probability)
- `TRANSMISSION_EFFICIENCY` (modifies a probability)
- `VERTICAL_TRANSMISSION` (modifies a probability)

### 28 February 2025
- Added missing `vertical_transmission_coefficient` to `PopulationType`, changed
vertical transmission names and parameter structure to conform with other
parameters
- Code for generating data demonstrating `flexle.jl` functionality (CLM)
- Added variables to control impact of responses on inoculum and vertical
transmission, but did not implement their use in the simulation: decided to
create a new kind of dictionary that contains coefficient functions that do not
correspond to events in the coefficient matrix but are still found in
`PathogenType` and `ResponseType`, this will reduce the number of repeated
functions we need to have in `immunity.jl` handling each one
- `FlexleSampler` weight update x10 performance improvements via: single
`FlexleSampler.sum` update; explicit + more specific type declarations (incl.
local variables); new `Dict` field in `FlexLevel` storing positions of items
in `indices` vector; and new `maxLevelWeight` function (CLM)

### 27 February 2025
- Added `model.time`
- Added parents to `Host`
- Added `birth_time` to `Host`, `Response`, and `Pathogen`
- Added option for Newick to show real time rather than Hamming distance
- Added `approxZero`, fixed bug in floating point error correction
- Fixed missing data error in `saveComposition`

I thought that individuals with immunity against a strain were not being
protected from coinfections with that strain, but that was only because we
weren't specifying an adequate `reactivityCoefficient` function.

### 26 February 2025
- Fixed bug where changes in contact receive weights were not correctly
propagating into contact rates (missing base population coefficient
multiplication) (PCR)
- Fixed floating point error that leads to very small but nonzero contact rate
(when it should be zero) (PCR)
- Fixed response acquisition bug (wrong constant) (PCR)
- Fixed response loss bug (wrong constant) (PCR)
- Fixed bugs in winner takes all immunity functions, which weren't taking the
`Response` with highest reactivity (PCR)
- Added `reactivityCoefficient`, `infectionCoefficient`, and
`responseStaticCoefficient` functions to streamline immunity and weight
functions (PCR)
- Set error tolerance constant to 1e-9 (PCR)
- Added error tolerance checking at model rate level (PCR)
- Moved `generateGamete` and `generateZygote` functions up in file (PCR)

### 25 February 2025
- Removed intrahost fitness coefficient from `PopulationType` since it was
never used (PCR)
- Changed default parameters to be all zeros for base coefficient/rates and
all maps to one for coefficient functions; Booleans were left as the most
common usage case (PCR)
- Reorganized examples into folders (PCR)
- Separated initializer functions for Type structs into a new file, `setup.jl`
- Renamed `parameters.jl` to `default.jl` (PCR)
- Fixed bug renormalizing the same random number with recombination events
during transmission (PCR)
- Changed structure of `developImmunity` arguments (PCR)
- Added logic to verify whether `Response` entities have actual `Pathogen`
entities associated, or just `nothing` (PCR)
- Added `deNovoResponse` as default `developImmunity` function (PCR)

KNOWN ISSUE: runtime didn't finish when I tried with
`response_acquisition_coefficient=1.0,`!!!

### 24 February 2025
- Added proportional fitness function for pathogen fractions (PCR)
- Debugged problems with complex weight calculation (contact, transition,
receive contact, receive transition), particularly (1) when adding hosts in
bulk, (2) when removing hosts, (3) in the main `hostWeights!` function (PCR)
- Write test function for `flexle.jl` user functions (CLM)
- Fix `flexle.jl` bugs mostly involving improper handling of zero-weight entries,
plus undeclared variable in `addToFlexleSampler!` as noted by PCR (CLM)

### 23 February 2025
- Remove function that plots compartments from output object to standardize
all plotting functions as working with `DataFrames` (PCR)
- Fix bugs in `compartmentPlot` function so that it works with `DataFrame`
input (PCR)
- Add total calculation to `saveComposition` (PCR)
- Make composition plots (PCR)
- Make Newick format phylogenies of pathogens, hosts, and responses (PCR)
- Plot phylogenies using the `NewickTree.jl` package--might be easier to plot
the newicks online though (PCR)
- Fix variable type bug in `zeroTruncatedPoisson` (PCR)
- Fix missing `transition!` function in event function constant vector (PCR)
- Change weight calculation at top level to minimize imprecision (overflow?)
errors (PCR)
- Removed `FunctionWrapper` type definitions in arguments of setup functions
because they can't handle passing functions as arguments apparently... hopefully
the type definitions at the struct level keep performance okay; the flamegraph
still doesn't show any non-compiled code at least (PCR)

### 22 February 2025
- Moved `logbounds` to `flexle.jl` (PCR)
- Type-defined `zeroTruncatedPoisson` and `catcol` (PCR)
- Fixed small bug in host weight calculation from adding host genomes (PCR)
- Separated out `mutateSequence` function (PCR)
- Add sexual reproduction option for `birth` (PCR)
- Added parameters for mutation and recombination of `Host` genome sequences
(PCR)
- Changed `mutateSequence` to account for chromosome separators (PCR)

Note: I think flexle's `addToFlexleSampler!` has an undeclared variable `i`

### 21 February 2025
- decided on `flexle.jl` API, functions written but not yet fully tested (CLM)
- Added genome sequence to `Host` (PCR)
- Made `Response` entities store the `Host` genome sequence they were derived in
(PCR)
- Made functions in `ResponseType` take host genome sequences into account (PCR)

### 20 February 2025
- `flexle.jl` additions including sampling; updating sampler `weights`; and
sample testing (CLM)
- Removed commented code for harmonic means and alternate algorithms in
immunity (PCR)
- Small bug fix in `developResponse` to pass correct arguments (PCR)
- Added a commented out `birth!` function to be used as a starting point for
sexual reproduction
- Separated out `recombineSequences` function (PCR)
- Changed `weightedResponse` and `infectionProbability` defaults to not be
the arithmetic means, but rather winner takes all style functions (PCR)

I decided that these `weightedResponse` and `infectionProbability` functions
must be imprinting-like functions. The simplest is winner takes all, but other
options could include top N respondents or some sort of a competition-based
algorithm. This could be something like this: you have a set amount of
resource for your total response (or response to a specific epitope,
corresponding to a `ResponseType`) and you allocate that resource to each
`Response` (of that `ResponseType`) starting from the one with the strongest
`weightedResponse` and going as far down the list as you still have resource
left to allocate. I don't love that specific idea though, not very mechanistic.
Regardless, as a default, it's best to do winner takes all, just like we're doing
with the default for intrahost pathogen fractions. For the "advanced" immunity
expansion/DLC, we might have different classes of `ResponseTypes` corresponding
to naive B cells, memory B cells, and long-lived plasmablasts, and at that
point we can rethink the best algorithm to use.

Based on intuition from this paper Amanda Elyssa Ruiz posted on Slack:
[10.1016/j.celrep.2025.115334](https://doi.org/10.1016/j.celrep.2025.115334)
which (copying from my Slack message) "vibes with the point Rachel \[Hecht]
raised during lab meeting that (if I'm remembering correctly) during
flu imprinting, when people look at GCs, they see a lot of B cell diversity,
but when they look at serum antibodies, they just see imprinted responses.
It vibes in the sense that this also suggests imprinting doesn't occur at
the level of selecting B cells that undergo maturation, but at the level of
selecting B cells to produce plasmablasts? If I understand things correctly?"

TODO:
- Add `Host` genome sequences, pass them to `developResponse`
- Add `Host` sexual reproduction, with second parent selection and
recombination during `Host` birth events

### 19 February 2025
- Finished composition dataframe (PCR)
- Bug fix on history dataframe saving (PCR)

### 18 February 2025
- Added data csv export of compartments data (PCR)
- Added csv plotting of compartments data (PCR)
- Added option of `nothing` for `Response.matured.pathogen` (PCR)
- Added data csv export of full data, can't believe it worked easily (PCR)
- Change `Host` sampling for history to a random sample of hosts fixed over
time (PCR)
- Started work on composition dataframe creation (PCR)

Next: composition dataframe and plot, then immunity

### 17 February 2025
- Small bug fix: `mean_recombination_crossovers` type corrected (PCR)
- Changed argument order in Type structs to put coefficients last within
their section (PCR)
- Finish model setup parameter initialization functions, default values (PCR)
- Extend `FlexleSampler` levels to include higher and lower bounds (CLM)
- `FlexleSampler` testing (CLM)
- Replace abstract `Function` type definitions in structs and function
arguments with `FunctionWrapper` functions with defined outputs and inputs (PCR)
- Type-define default argument functions (PCR)
- Identified performance issue using flamegraph in which mutable static vector
defining `inocula` within `hostContact!` resulted in runtime function resolution;
made `inocula` into a regular `Vector` and this removed the runtime resolution
of the calls (but addded orange bars on the flamegraph in what I'd consider
unrelated parts of the code, the return lines of `choosePathogen!` and
`chooseHost!`--worth keeping an eye on, maybe with `inocula` of length > 1 we
will be able to spot performance problems and `MVector` will be better) (PCR)

At this point, the flamegraph shows no runtime (red) calls!

### 16 February 2025
- Bug fix in compartment variable handling (PCR)
- New `Output` struct for simulation output (PCR)
- Population compartment plot (PCR)
- Moved sandbox file to examples folder (PCR)
- Default model setup functions (PCR)
- Add/remove/move indices to/from/between `FlexleSampler` levels (CLM)

Next:
- Finish model setup functions
- Immunity draft

### 15 February 2025
- Moved FlexLev/Flexle code into a new file for modularization of the technique
(PCR)
- Added `StaticHost` for more lightweight storage of host info in simulation
history (PCR)
- Added history tracking in two ways: accurate counts of compartment variables
(uninfected naive, infected naive, uninfected immune, infected immune, and dead
hosts) in each population at each indicated time point, and a sample of hosts
(each stored as `StaticHost`) from each population at each time point (PCR)
- Reworked `Flexle` sampling to use arrays as described below (CLM)

### 14 February 2025
- Added references to actual parent entities in `Pathogen` and `Response` using
CLM's `Union{X, Nothing}` trick (PCR)
- Implemented flexible level + rejection sampling (current working titles "FlexLev
sampling", or perhaps "Flexle sampling") as described in PCR's note yesterday.
However, in discussions today re: the relative benefits of using a linked list
structure for `FlexLevel` entities in a `FlexlevSampler` as opposed to an array,
we decided that the primary LL advantage—namely the ability to easily add new
`FlexLevel` without copying data—is likely to be outweighed by both (1) that arrays
allow for constant-time indexing, meaning adding, removing, and moving weights
to/from/between levels (actions which are all likely to be far more common than
creating new levels) are much faster for arrays than for linked lists, and (2) that
the "add level" operation is actually still _O(n)_ in the number of levels for linked
lists when levels are added in the middle of the list, thereby partially negating the
linked list's primary benefit. (As a side note, arrays actually lose their constant-time
indexing property for our use case if you allow for non-contiguous `FlexLevel`s in a
`FlexlevSampler`. However, we have decided for now at least that the additional
space/time to enforce level contiguousness is an appropriate price to pay for
constant-time indexing.) All of that is to say that
*flexle sampling as currently implemented using linked lists will not remain*.
CLM will rework the implementation to use arrays instead. Nevertheless, this iteration
of flexle sampling is being committed for posterity. (CLM)

### 13 February 2025
No changes yet, but based on (unnecessarily, sorry) lengthy discussion, we conclude
that the optimal implementation of random weighted sampling of `Host` entities is
a modified version of the binary level rejection sampling technique described by
[Aaron DeFazio](https://www.aarondefazio.com/tangentially/?p=58) based on work
from Donald Knuth and others ("The technique used is not novel, indeed it is
based on publications from the 1960s" says DeFazio--but who?!).

We will implement rejection sampling with binary levels as described by DeFazio &
1960s, including optimizations such as random number renormalizing as done in Opqua
already and also explained by DeFazio in the "Enhancements" section. We discussed
whether we should include in the algorithm a check on whether to add a new level
at the top or bottom of the vector of levels, in case the ratio between the largest
and smallest weights exceeded the 2^20-fold recommended by DeFazio (who argues for
20 intervals, given that the dynamic range of weight values is usually less than
1 million, but nothing guarantees this). We decided to do this, and we might as
well also remove unnecessary levels at the top or bottom of the level list if they
become empty of weights. In order to most efficiently add and remove levels both at
the start and end of the list of levels, CLM suggested we use a linked list. We
also discussed whether rejection sampling of the list of levels itself would
outperform linear searching our way through it, as implemented by DeFazio, but we're
unsure, given the potentially lopsided differences in the sums of weights between
different levels. We decided we will stick with linear search across levels.

We are pretty confident this is a good idea for `Host` sampling given the number of
hosts in a `Population`. We suspect the current linear search algorithm implemented in
`randChoose` is optimal for sampling event types, given there are only 9 events:
according to tests done by CLM in Julia, a uniformly-balanced weights vector is
sampled faster by linear search than by simple rejection sampling if there are fewer
than 10 weights being sampled* (PCR is now unsure of this, see note below). We are
unsure whether sampling of `Pathogen` and `Response` entities within a `Host` is best
done with linear search or rejection sampling given that the distribution of the
number of pathogens and responses is both variable and skewed in favor of smaller
numbers. Because of the latter, we will remain on linear search for the time being.

*Bear in mind that sampling from a uniform weight array represents the best possible
scenario for simple rejection sampling, since the algorithm reduces to indexing into
a particular position in the vector based on a random integer and then always
accepting the given element regardless of the second random number that accepts or
rejects the chosen weight/element. However, I now remember that this second random
number can be obtained by renormalization instead of fresh regeneration. CLM did not
implement this detail in the test, so maybe the vector size at which rejection
sampling becomes faster than linear search is closer to 5 than 10? This makes me less
confident that linear search is preferrable now; perhaps testing binary level
rejection sampling for events and even pathogens and responses is actually worth it.

### 12 February 2025
- Removed `total_hosts`, replaced with instances of `length(population.hosts)` (PCR)
- Optimized `Host` initialization at startup (CLM)
- Made `event_functions` into constant `EVENT_FUNCTIONS` (PCR)

### 11 February 2025
- Added interventions (PCR)
- Added flamegraphs to sandbox file (PCR)
- Defined a couple of static array types within `events.jl` (PCR)
- Made new matrices in `Population` with sums of host weights pre-multiplied by
coefficients, to be used in `chooseHost` replacing the inline, broadcasted,
element-wise multiplication by the base coefficient; resulted in a ~4X speedup and
garbage collection time going from ~18% to ~1.5%, even with the speedup--the problem
was found with the dual wonders of flame graphs and CLM (PCR)
- Defined type of `AbstractArray` in `randChoose` (PCR)

There was a bug in the use of parenthesis for `@views` in `choice.jl`, very
minor effect on performance. The fix was made obsolete by removing the broadcasted
multiplication anyway.

Next is saving history--a clean implementation of this would take advantage of
the intervention framework to save whatever is needed whenever needed, but might
result in performance issues given the type instability of interventions, so maybe
a standalone implementation of sampling and saving might be best?

### 9 February 2025
- Small bug fix in `transition!` (PCR)

### 5 February 2025
- Added pathogen fraction to calculation of vertical transmission probability (PCR)
- Moved code to add a `Host` to a `Population` into a separate function (PCR)
- Optimizations including: add types to global variables; match literal types to
their context; extend use of `@views`; cache emitting `Host` in `hostContact!` (CLM)
- Backtrack 2 optimizations (one `@views` use with no slicing, one comparison between
ints that did not require using a float zero) (PCR)
- Added `removeHostFromPopulation` and death event (PCR)
- Added transition event (PCR)

All events implemented!!! but not tested lol

TODO:
- function for running and benchmarking repeated simulations; current setup
has a global `model` variable which I (CLM) believe retains its state after
simulation _n_ at the beginning of simulation _n+1_, affecting both simulation
trajectory and runtime
- Test and debug the **** out of everything (all parameters and events,
particularly events other than clearance and contact)

Dev roadmap:
- implement interventions
- simulation history saving
- input parameters
- output data + plotting
- immunity DLC

### 4 February 2025
- Added `constant_contact_density` as a parameter of `PopulationType` to specify
whether or not to normalize by host population size in calculation of
receive contact weights, equivalently, added `constant_transition_density`
for receive transition weights (`constant_transition_density` should by
default be `true`; can't imagine otherwise)
- Add response loss event (not debugged)
- Add birth event, along with response inheritance and vertical transmission
parameters (not debugged)

### 3 February 2025
Significant overhaul of weights calculations in order to accommodate correct
inter- and self-population contact

More than just this, but some notes I took:
- Change weights calculation for contacts at population level so that it
incorporates infected and total host numbers as well as inter-population
contacts
- Change weight propagation so that event weights are not propagated above
population level if the number of hosts is zero, and if they are zero, the
corresponding population weight is set to zero and propagated
- Implement contacts as inter- or intra-population contacts

### 31 January 2025
- Moved `hostContact!` down to section with events
- Renamed `MIGRATION` to `TRANSITION` and `migration_fractions` to
`transition_rates`
- Changed name of `PopulationParameters` to `PopulationType`
- `pathogenFractions!` is now `pathogenFractionsWinnerTakesAll`,
a parameter of `PopulationType`; default function is stored in
new file, `intrahost.jl`
- Made `weightedResponse` and `infectionProbability` into
`weightedResponseArithmeticMean` and `infectionProbabilityArithmeticMean`,
parameters of `PopulationType`; default functions are stored in
`immunity.jl`
- Moved `immunity.jl` and `intrahost.jl` into their corresponding
directories within the new `extensions` directory
- Add inoculum, mutation, recombination population-level coefficient
parameters
- Rename `population_contact_sum` to `contact_sum`

TODO: Implement inter-population contacts such that self-contact is an option
and we just have a single contacts function for all cases; receive contact
weights should be handled like contact weights, one per possible population
contact relationship


### 29 January 2025
- Removed `Class` and everything associated to it
- Made a single `CONTACT` event

### 23 January 2025
- Small bug fix in `newHost!`
- Track `intra_population_contact_sum` within `Population`
- Made intra-population contact rate at the `Population` level account for
receive contact and total host populations
- Changes in event modifier weights are now computed before events so that
receive weights are computed prior to contact rates
- Make changes in receive weights propagate into contact rates at the
`Population` level (different from what I described yesterday, but I like
this better)
- Fix weight propagation

#TODO: continue debugging contacts, not convinced simulation trajectories are
reasonable... actually, I think it works
If we can make inter-population contact rates depend on receive rates the way
intra-population contact works, we should consider two major changes in
architecture: (1) getting rid of `Classes` completely since `Populations` can
model anything a `Class` does, and (2) getting rid of the distinction between
intra- and inter-population contacts entirely, since an inter-population contact
with the same population achieves the same objective.

### 22 January 2025
- Make `Class`, `Population`, and `Model` receive weight matrices include weights
for all coefficients except for fitness
- Change `intraPopulationContact!` such that host 2 is sampled using receive
contact coefficients (instead of intrapopulation contact coefficients)
- Minor debug in `hostWeightsReceive!`, missing a couple lines
- Created `attemptInfection` function to handle checking for infection likelihod,
etc. before successfully adding a `Pathogen` to a `Host` so that it can be
called by `intraPopulationContact`, `interPopulationContact`, `establishMutant`,
`establishRecombinant`, and `birth` (if there is vertical transmission)
- Move `RECOMBINATION_PER_REPLICATION` and `MUTATION_PER_REPLICATION`
to named fields of `PathogenType`, add corresponding coefficient functions (and a
function for effective inoculum), make corresponding fields within `Pathogen`;
don't bother with the specific or static coefficient functions  in `Response`,
they don't matter; considered moving fitness to a separate field too but because
fitness is actually used in sampling, this is incorrect (rationale is the other
two coefficients are not used to sample entities)
- Changed computation of mutation and recombination probability during
transmission to be based on Poisson process means
- Moved contact algorithm into its own function to be called by
`intraPopulationContact` and `interPopulationContact` with a parameter specifying
whether or not to resample populations (#TODO: sampling of population 2 needs to be
done according to contact rates from population 1)
- Added propagation functions that go from the population level to the event level

KNOWN ISSUE: contact rate computation is incorrect; we have to multiply the sum of
contact rates by the sum of receive weights inside the `Class` level and divide by
the total number of hosts inside the `Population` level. This is a problem for the
current propagation architecture because it can handle changes in a single
variable at a time--to account for a change in contact weight and a change in
receive weight stemming from a single event, it currently has to go through
propagation cascades twice. To solve this, we will add a series of propagation
functions specifically for contact rates that take not only the change in contact
rate as a parameter, but also the simultaneous change in receive rate.

- Relatedly, we added `total_hosts` and `receive_contact_sum` fields to
`Population` to be used in the computation of contact rate and modified
`newHost!` to change both `total_hosts` and contact weights at the `Population`
level correspondingly

### 21 January 2025
- Some simple syntax debugging in events functions
- Important debugging in hostWeights computation, some weight matrix indexes were
missing and weight changes for propagation were not computed correctly for pathogen
and response events
- Make populations retain `RECEIVE_CONTACT` weights in `class_weights`
- Added `immunity.jl` to contain functions to simulate immune responses using
`Response` elements, added `immunityProbability` as a function that computes
probability that a `Host` is immune to a specific infecting `Pathogen`
- Added `immunityProbability` field to `ResponseType` that computes the
probability that a `Host` is immune to a specific infecting `Pathogen` thanks to
a `Response` of this type
- Changed `intraPopulationContact!` to check for immunity with `immunityProbability`
before adding a `Pathogen`
- Fix weight/rate computation to modify rate sum (wasn't happening)
- Replaced `sample()` with `chooseRandom()` in `intraPopulationContact`, Connor
benchmarked and said it was ~5X faster, at least sometimes?
- Renamed `immunityProbability` field within `ResponseType` as `infectionCoefficient`,
changed mechanism such that a coefficient of `1.0` is complete susceptibility and
`0.0` is complete immunity
- Added `reactivityCoefficient` as a field of `ResponseType` to handle how much a
`Response` contributes to a reaction against a specific pathogen
- Change sums of `Response` specific coefficients to weighted arithmetic means,
added `reactivityCoefficient` to calculation as weight, created a `weightedResponse`
function to compute all specific weighted mean responses (`infectionCoefficient`
weighted mean is computed separately); thought about whether these should be
weighted arithmetic or harmonic means, I think arithmetic (the work is being done
by the `reactivityCoefficient` function defining weights)
- Added logic to make recombination weights/probability zero when less than two
strains

Considered renaming `pathogen_fractions` to `pathogen_populations` and
`pathogenFractions` to `pathogenPopulations`, change the latter's behavior such
that the top pathogen's population is reported as its coefficient product, not just
`1.0`. Decided for sure against it; `pathogen_populations` captures competition,
actual pathogen-specific strength/influence on events is captured by other functions.

KNOWN ISSUES: infection events still don't seem to result in infection. This is a
problem for tomorrow. First time simulation runs though!

### 20 January 2025
- Rename `acquireResponse` to `developResponse` to avoid ambiguity about when a
`Response` is actually added to a `Host`
- Added a quasi zero-truncated Poisson distribution function as used in Python Opqua;
worth looking into whether an actual zero-truncated Poisson is preferrable
- Added `num_recombination_crossovers` and `mean_effective_inoculum` parameters to
`PathogenType`
(as an aside, I have decided to use Poisson distributions to model inoculum size
rather than negative binomials for simplicity and due to the fact that the increased
variability in mean inoculum size seen in Sobel Leonard _et al._ (2017) may have been
inflated due to the technical contaminations in Poon _et al._ (2016),
(https://doi.org/10.1128/jvi.00936-19)[https://doi.org/10.1128/jvi.00936-19])
- Change `pathogens` dictionary in `Class` to use genomes as keys and
correspondingly changed `responses` dictionary to use tuples of genomes as keys;
removed `id` field from both `Pathogen` and `Response` and `pathogens_idx` from
`Class`
- Rename `RECOMBINATION` event as `RECOMBINANT_ESTABLISHMENT`, rename
`MUTATION_AT_CONTACT` as `MUTATION_PER_REPLICATION`, add a new event choice
coefficient `RECOMBINATION_PER_REPLICATION`
- Added `binomial()` to generate binomial random numbers, faster than
`Distributions.jl`
- Added event functions for mutant establishment, pathogen clearance, response
acquisition, recombinant establishment, and intra-population contact--debugging
needed for these

### 19 January 2025
- Changed `pathogens` and `responses` vectors in `Class` to be dictionaries instead,
with keys corresponding to the assigned `Pathogen` integer code in the case of
`pathogens` and a tuple with the codes for imprinted and matured pathogens
in the case of `responses`
- Add `parent` as a new field of `Response`, storing a tuple referencing the
`Response` that this `Response` was derived from within its `Class`; this is only
useful for response lineage tracing, but not the simulation itself (at least at
the moment), so we may remove it
- Added base function to create a mutant `Pathogen`
- Added vector to track population receive weight sums in `Model`
- Changed weight calculation to update sum of weights in next level of hierarchy as
you iterate through the current level
- Added `Choice.jl` with functions to sample different entities based on their weights

### 18 January 2025
- Changed language such that all concepts relating to "immunity" now use "response"
- Renamed functions to add and remove pathogens and immunities from hosts
- Added new event, `RESPONSE_ACQUISITION`
- Removed `immunodominance` from `ResponseType`; this can be handled through
user-defined functions specifying how `Response` elements are acquired during
an `RESPONSE_ACQUISITION` event
- Added `PathogenType` to store `Pathogen`-specific parameters

### 17 Jan 2025
- Change `type` field within `Immunity` to be a reference to the `ImmunityType`,
not an index

### 12 Jan 2025
- Fixed many weight computation bugs

### 8 Jan 2025
- Added infection, immunization, clearance, deimmunization functions
- Changed all `Pathogen`, `Immunity`, and `ImmuneType` arrays to be arrays of
(references to) the entities, rather than arrays of their integer indexes
TODO: revise coefficient/weight/event nomenclature, clean up `weights.jl` in particular

### 5 Jan 2025
- Debugged `Host` weight initialization
TODO: Adding infection and immunization functions

### 4 Jan 2025
- Reorganized functions into new file structure
- Add and first debug struct creation methods
- Change `Pathogen` and `Immunity` entity storage within `Class` to use a simple vector
of entities, rather than a dictionary with their associated ID
- Change `Pathogen` and `Immunity` entity storage within `Host` to use vectors of entities
(references to them), rather than their indexes in another matrix
To do: Code to sample events and entities based on rates and weights

### 31 Dec 2024
Central change here is implementing a system for handling immunity through
`Immunity` entities that store references to the `Pathogen` sequences that were
imprinted and affinity matured to generate them and a reference to an `ImmunityType`
entity that defines how the immunity affects coefficients based on its imprinted
and matured sequences, as well as any infecting pathogens.
More changes:
- Moved struct definitions into `structs.jl` to compile before methods
- Changed `Pathogen` and `Immunity` vectors within `Hosts` to be vectors of keys
 used to access those objects within `Class` dictionaries instead
- Remove immune priority and the previous scheme used for weighting the effect of
different immunities; now, effect of immunities is not dependent on other pathogens
present
- Added `ImmunityType` struct to store "epitope types"; `Immunity` entities now
 have an associated `ImmunityType` that defines their coefficient functions and
 immunodominance; `ImmunityType` entities are stored inside `ClassParameters`
- Added `acquireImmunities` function field to `ClassParameters` to define which and
 how many `Immunity` entities (with corresponding `ImmunityTypes`) are acquired during
 immunization
- Changed `Immunity` to store the index of a `Pathogen` within the pathogen vector
in `Class` rather than storing a sequence
- Added coefficient vector to `Immunity`, make those coefficients affect weights
(these are coefficients that are "statically" active, i.e. they don't compare
sequences)
- Changed `Immunity` to store both imprinted and matured `Pathogen` indexes, then
change immunity functions such that they use up to three pathogen sequences (imprinted,
matured, and infecting) instead of only comparing two (immunized and infecting); in
this logic, imprinted immunity would have an imprinted pathogen sequence, but not a
matured one


### 16 Dec 2024
- Some quick bug fixes
- Added quick rundown of model structure to README
To do: `Immunity` is currently just a sequence, can be just a reference to the corresponding 
pathogen.

### 15 Dec 2024
- Made fitness into a pathogen-level coefficient
- Renamed "transmission" language to "contact"
- Distinguished between `INTRA_POPULATION_CONTACT` and `INTER_POPULATION_CONTACT`
- Changed transmission-based mutation into an event modifier coefficient rather
  than its own independent event

First draft of the following three:
- Added `Class` matrices for sampling pathogen, immunity, and host event weights according to event
- Added `Population` matrix for sampling class event weights according to event
- Added `Model` matrix for sampling population event weights according to event

Next: Code `Model` vector for sampling total event rates with correct rate computation,
add code for sampling class change and migration destination using corresponding rates in class.

### 13 Dec 2024
First draft of `Pathogen`, `Immunity`, and `Host` coefficients, fractions, and weights.
Next up are `Class` and real rates.

### 10 Dec 2024
First attempt at structs and event vector structure.

### 9 Dec 2024
Created directory and file structure.
