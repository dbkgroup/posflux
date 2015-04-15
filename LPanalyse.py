import math

def _extract_flux_patterns_from_PositiveFBA(LPsol,num_reactions):
	varValues = LPsol.getVars()
	flux = []
	for tic in range(4*num_reactions,5*num_reactions):
		component = varValues[tic]
		flux.append(component.X)
	return(flux)
	
def _extract_flux_patterns_from_SuperDaaave(LPsol,num_reactions):
	#Equivalent to fetching variables e
	fluxes = LPsol.getVars()
	fluxVector = []
	for tic in range(num_reactions):
		fluxObj = fluxes[tic+(4*num_reactions)]
		#Tidy up some weird sign errors from Gurobi (minus zero?)
		f = fluxObj.X
		f_squared = math.pow(f, 2)
		f_sqrt = round(math.sqrt(f_squared),4)
		fluxVector.append(f_sqrt)
		
	return(fluxVector)

def _extract_flux_patterns_from_ComparisonDaaave(LPsol,num_reactions):
	varValues = LPsol.getVars()
	fluxVector1 = []
	for tic in range(4*num_reactions, 5*num_reactions):
		fluxObj = varValues[tic]
		#Tidy up some weird sign errors from Gurobi (minus zero?)
		f = fluxObj.X
		f_squared = math.pow(f, 2)
		f_sqrt = round(math.sqrt(f_squared),4)
		fluxVector1.append(f_sqrt)

	fluxVector2 = []
	for tic in range(9*num_reactions, 10*num_reactions):
		fluxObj = varValues[tic]
		#Tidy up some weird sign errors from Gurobi (minus zero?)
		f = fluxObj.X
		f_squared = math.pow(f, 2)
		f_sqrt = round(math.sqrt(f_squared),4)
		fluxVector2.append(f_sqrt)
		
	return(fluxVector1,fluxVector2)