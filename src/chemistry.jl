"""
O2+ dissociative recombination chemistry.

The coefficient and four fixed branches follow the project configuration.
Reaction rates are converted to SI units before particle sampling.
"""

"""O2+ dissociative recombination coefficient in m^3 s^-1."""
function dissociative_recombination_coefficient(Te_K::Real)
    Te_K > 0 || throw(DomainError(Te_K, "Electron temperature must be positive"))
    coefficient_cm3_s = Te_K <= 1200 ?
        1.95e-7 * (300 / Te_K)^0.70 :
        7.39e-8 * (1200 / Te_K)^0.56
    coefficient_cm3_s * 1e-6
end

"""Production rate of both hot O products, in m^-3 s^-1."""
hot_o_production_rate(ne_m3::Real, nO2p_m3::Real, Te_K::Real) =
    2ne_m3 * nO2p_m3 * dissociative_recombination_coefficient(Te_K)

function load_reaction_branches(path::AbstractString)
    input = TOML.parsefile(path)
    branches = [DRBranch(
        item["products"], item["release_energy_eV"], item["probability"],
    ) for item in input["branch"]]
    isapprox(sum(x.probability for x in branches), 1.0; atol=1e-12) ||
        error("Dissociative recombination branching probabilities must sum to one")
    branches
end
