import SBMLRead as sbml_read
import Positive_Constraint_based_Modelling as pcbm
import LPSolver as solver
import LPanalyse as analyse

model = 'data/yeast_5.21_MCISB.xml' 

Smatrix, num_species, num_reactions, reaction_names, lbs, ubs, gene_associations = sbml_read.retrieve_information_from_SBML(model)

#Force glucose into the system (note the exchange reaction is written in the above model as glucose moving from inside the cell to outside, so the flux is negative in the SBML).
glucoseIndex = reaction_names.index('D-glucose exchange')
lbs[glucoseIndex] = -16.5
ubs[glucoseIndex] = -16.5

#Set the the growth reaction as the objective
growthIndex = reaction_names.index('growth')
objective = {}
objective[growthIndex] = -1


LHS,equalities,RHS,variable_types,lbs,ubs,objective = pcbm._buildPositiveFBA(Smatrix,lbs,ubs,objective)									#Build matrix problem

LP = solver._create_GurobiLP(LHS,equalities,RHS,variable_types,lbs,ubs,objective)														#Build solver problem

LPsolved = solver._run_GurobiLP(LP)																										#Solve LP

flux = analyse._extract_flux_patterns_from_PositiveFBA(LPsolved,num_reactions)															#Extract flux pattern

for i in range(len(flux)):
	print round(flux[i],3), '\t', reaction_names[i]