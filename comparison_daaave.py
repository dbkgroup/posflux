'''
(c) University of Liverpool 2019

All rights reserved.

The ComparisonDaaave algorithm extends the principle of the original Daaaaave
algorithm that seeks flux patterns which correlate with measured
transcript/protein abundances.
Instead ComparisonDaaave algorithm seeks two flux patterns, the ratios of which
correlate with measured ratios of proteins (as from iTRAQ, for example).

The algorithm requires a model in SBML, experimental data in TSV format
(gene\tmean) - note stdevs are not yet accounted for within the objective,
and...

...the model must be forced into life by pushing flux through reactions.
The formulation of the problem is such that the 'best' solution is to select
flux1 as zero as this allows flux2 to also be zero and the ratio is ignored.
'''
import lp_analyse
import lp_solver
import matplotlib.pyplot as plt
import positive_constraint_based_modelling
import sbml_data_mapper


def main():
    '''main method.'''
    model = 'data/yeast_5.21_MCISB.xml'
    data = 'data/yeastComparison.txt'

    print('- Fetch model and map in relative expression data')
    s_matrix, _, num_reactions, reaction_names, lbs, ubs, \
        rxn_mean_ratio \
        = sbml_data_mapper._read_and_map_sbml_and_relative_data(model, data)

    # Force flux through the system
    lbs1 = lbs
    ubs1 = ubs
    lbs2 = lbs
    ubs2 = ubs

    glucose_index = reaction_names.index('D-glucose exchange')
    glucose_flux_value1 = -16.5
    glucose_flux_value2 = -5.00
    lbs1[glucose_index] = glucose_flux_value1
    ubs1[glucose_index] = glucose_flux_value1
    lbs2[glucose_index] = glucose_flux_value2
    ubs2[glucose_index] = glucose_flux_value2

    growth_index = reaction_names.index('growth')
    growth_flux_value1 = 1.00
    growth_flux_value2 = 0.50
    lbs1[growth_index] = growth_flux_value1
    ubs1[growth_index] = growth_flux_value1
    lbs2[growth_index] = growth_flux_value2
    ubs2[growth_index] = growth_flux_value2

    # Build ComparisonDaaave
    print('- Build ComparisonDaaave problem')

    # Build matrix problem
    lhs, equalities, rhs, variable_types, lbs, ubs, objective = \
        positive_constraint_based_modelling.build_comparison_daaave(
            s_matrix, lbs1, ubs1, lbs2, ubs2, rxn_mean_ratio)

    print('- Create Gurobi instance of problem')

    # Build solver problem
    lp = lp_solver._create_GurobiLP(
        lhs, equalities, rhs, variable_types, lbs, ubs, objective)

    print('- Solve ComparisonDaaave')
    LPsolved = lp_solver._run_GurobiLP(lp)  # Solve LP

    flux1, flux2 = lp_analyse._extract_flux_patterns_from_ComparisonDaaave(
        LPsolved, num_reactions)  # Extract flux patterns for models 1 and 2

    x = []
    y = []

    for i in range(len(rxn_mean_ratio)):
        if rxn_mean_ratio[i] > -1:
            if flux1[i] != 0 and flux2[i] != 0:
                f1 = flux1[i]
                f2 = flux2[i]
                fluxRatio = f2 / f1
                x.append(rxn_mean_ratio[i])
                y.append(fluxRatio)
                # x.append(math.log10(rxn_mean_ratio[i]))
                # y.append(math.log10(fluxRatio))

    # Plotting some results
    plt.scatter(x, y)
    plt.title('Protein Ratio vs Nearest Possible Flux Ratio')
    plt.xlabel('log(Protein Ratio)')
    plt.ylabel('log(Flux Ratio)')
    plt.show()
