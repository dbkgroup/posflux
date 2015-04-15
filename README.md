# Positive Flux Modelling
Constraint-based metabolic modelling with only positive fluxes.

Python coding of algorithms for constraint-based modelling of metabolism with reactions re-written so that only positive fluxes are considered. These positive flux models consist of +S (the stoichiometric matrix), -S (all the reactions written in reverse), some fiddling with the reaction bounds to preserve original reaction directions, and sets of binary variables to ensure only one version of each reaction (+S or -S version) is active at any one time. Technically they are all MILPs, not simple linear programs.

NOTE: Gurobi is used to solve all LPs. A standard LP formulation is provided, so it should be relatively easy to extend the code to accommodate other solvers.

## Positive FBA

This will be slower than FBA implementations that do not worry about negative flux, but as it is at the core of the other algorithms it is being made available anyway.


## SuperDaaave

A high-speed version of the Daaaaave algorithm, written as one MILP rather than iterative LPs.
Daaaaave is a constraint-based algorithms that seeks flux patterns which correlate with expression data.

The original Daaaaave algorithm is described in:
Lee D, Smallbone K, Dunn WB, Murabito E, Winder CL, Kell DB, Mendes P, Swainston N (2012) "Improving metabolic flux predictions using absolute gene expression data" BMC Syst Biol 6:73
http://dx.doi.org/10.1186/1752-0509-6-73
  
The original Daaaaave algorithm is available at http://github.com/u003f/daaaaave/releases/tag/original

## ComparisonDaaaave (BETA)

Like SuperDaaave but for relative abundance data. Maximises correlation between relative abundance and relative fluxes.
