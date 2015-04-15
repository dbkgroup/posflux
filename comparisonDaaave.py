import SBMLDataMapper as sbml_data_mapper
import Positive_Constraint_based_Modelling as pcbm
import LPSolver as solver
import LPanalyse as analyse
import math
import matplotlib.pyplot as plt

'''
The ComparisonDaaave algorithm extends the principle of the original Daaaaave algorithm that seeks flux patterns which correlate with measured transcript/protein abundances.
Instead ComparisonDaaave algorithm seeks two flux patterns, the ratios of which correlate with measured ratios of proteins (as from iTRAQ, for example).

The algorithm requires a model in SBML, experimental data in TSV format (gene\tmean) - note stdevs are not yet accounted for within the objective, and...

...the model must be forced into life by pushing flux through reactions. The formulation of the problem is such that the 'best' solution is to select flux1 as zero as this allows flux2 to also be zero and the ratio is ignored.

'''

model = 'data/yeast_5.21_MCISB.xml' 
data = 'data/yeastComparison.txt'

print '- Fetch model and map in relative expression data'
Smatrix, num_species, num_reactions, reaction_names, lbs, ubs, rxn_mean_ratio = sbml_data_mapper._read_and_map_sbml_and_relative_data(model,data)

#Force flux through the system
lbs1 = lbs
ubs1 = ubs
lbs2 = lbs
ubs2 = ubs

glucoseIndex = reaction_names.index('D-glucose exchange')
glucoseFluxValue1 = -16.5
glucoseFluxValue2 = -5.00
lbs1[glucoseIndex] = glucoseFluxValue1
ubs1[glucoseIndex] = glucoseFluxValue1
lbs2[glucoseIndex] = glucoseFluxValue2
ubs2[glucoseIndex] = glucoseFluxValue2

growthIndex = reaction_names.index('growth')
growthFluxValue1 = 1.00
growthFluxValue2 = 0.50
lbs1[growthIndex] = growthFluxValue1
ubs1[growthIndex] = growthFluxValue1
lbs2[growthIndex] = growthFluxValue2
ubs2[growthIndex] = growthFluxValue2

#Build ComparisonDaaave
print '- Build ComparisonDaaave problem'
LHS,equalities,RHS,variable_types,lbs,ubs,objective = pcbm._buildComparisonDaaave(Smatrix,lbs1,ubs1,lbs2,ubs2,rxn_mean_ratio)			#Build matrix problem

print '- Create Gurobi instance of problem'
LP = solver._create_GurobiLP(LHS,equalities,RHS,variable_types,lbs,ubs,objective)														#Build solver problem

print '- Solve ComparisonDaaave'
LPsolved = solver._run_GurobiLP(LP)																										#Solve LP

flux1, flux2 = analyse._extract_flux_patterns_from_ComparisonDaaave(LPsolved,num_reactions)												#Extract flux patterns for models 1 and 2
x = []
y = []
for i in range(len(rxn_mean_ratio)):
	if rxn_mean_ratio[i] > -1:
			if flux1[i] != 0 and flux2[i] != 0:
				f1 = flux1[i]
				f2 = flux2[i]
				fluxRatio = f2/f1
				x.append(rxn_mean_ratio[i])
				y.append(fluxRatio)
				#x.append(math.log10(rxn_mean_ratio[i]))
				#y.append(math.log10(fluxRatio))

#Plotting some results
plt.scatter(x,y)
plt.title('Protein Ratio vs Nearest Possible Flux Ratio')
plt.xlabel('log(Protein Ratio)')
plt.ylabel('log(Flux Ratio)')
plt.show()
