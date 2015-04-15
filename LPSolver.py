import numpy as np
from scipy.sparse import csr_matrix, hstack, vstack
from scipy import sparse
import math
from gurobipy import *

def _run_GurobiLP(lp):
	lp.optimize()
	return(lp)

def _create_GurobiLP(LHS,equalities,RHS,variable_types,lbs,ubs,objective):
	'''
	LHS				is a CSR matrix of size c * v (constraints * variables) with nonzero coefficients
	equalities		is a list of size c containing equalities, inequalities
	RHS				is a list of size c of values (for now no functions are included)
	variable_types	is a list of variables types
	lbs,ubs			are lists of size v containing lower and upper bounds
	objective		is the dictionary of non-zero elements and their weights (virtually of size v, though in practice often smaller)
	'''
	#Create model
	lp = Model("LP")
	
	#Create variables, note all variables here are continuous
	print '\tCreate Gurobi variables'
	for varIndex in range(len(variable_types)):
		lb = lbs[varIndex]
		ub = ubs[varIndex]
		
		#Identify variable type
		if variable_types[varIndex] == 'continuous':
			var_type = GRB.CONTINUOUS
		elif variable_types[varIndex] == 'binary':
			var_type = GRB.BINARY
		else:
			var_type = 'unknown'
		
		#Identify whether objective
		if varIndex in objective:
			objective_value = objective.get(varIndex)
		else:
			objective_value = 0.0

		#Create variable
		lp.addVar(lb,ub,objective_value,var_type,'v'+str(varIndex))

	#Integrate new variables
	lp.update()

	lpvars = lp.getVars() #List of variable objects
	#Create constraints
	print '\tCreate Gurobi constraints'
	for consIndex in range(len(equalities)):
		#sys.stdout.write('\t\t' + str(consIndex) + '\r')
		#sys.stdout.flush()
		lin = []	#Coefficients
		refs = []	#Variables
		row = LHS.getrow(consIndex)
		_,vars = row.nonzero()
		vars = vars.tolist()
		lin = row.data
		for r in vars:
			coeff = row.getcol(r).toarray()[0][0]
			#lin.append(coeff)
			refs.append(lpvars[r])

		if equalities[consIndex] == '>':
			equality = GRB.GREATER_EQUAL
		elif equalities[consIndex] == '<':
			equality = GRB.LESS_EQUAL
		elif equalities[consIndex] == '=':
			equality = GRB.EQUAL
		else:
			equality = ''

		#print consIndex, lin, refs, equality, RHS[consIndex]
		lp.addConstr(LinExpr(lin,refs), equality, RHS[consIndex],'c'+str(consIndex))
	
	#Integrate new constraints
	lp.update()
	
	return(lp)