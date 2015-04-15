from SBMLContentHandler import SBMLContentHandler
import xml.sax
from scipy.sparse import dok_matrix, csr_matrix
import numpy as np

'''
Fetches information from SBML by;

Smatrix, num_species, num_reactions, reaction_names, lbs, ubs, gene_associations = retrieve_information_from_SBML_initialise_model(path_to_sbml_file,default_bound)

e.g.
Smatrix, num_species, num_reactions, reaction_names, lbs, ubs, gene_associations = retrieve_information_from_SBML_initialise_model('../Dave/mcisb-constraint-based-modelling/data/yeast/yeast_5.21_MCISB.xml',100000000)

Creates:
	S matrix = a stoichiometric matrix in sparse CSR format of num_species rows and num_reactions columns
	num_species = the number of species in the SBML
	num_reactions = the number of reactions in the SBML
	reaction_names = a list of length num_reactions containing the names of each reaction in the SBML
	lbs, ubs = lists of length num_reactions containing the lower and upper bounds for each reaction
	gene_associations = a list of length num_reactions containing the genes associated with each reaction
'''

def _initialise_model(sbml,default_bound=1000000):
    sbml_content_handler = SBMLContentHandler(default_bound)
    xml.sax.parse(open(sbml), sbml_content_handler)
    species_reactions, species_stoichiometries, num_species, num_reactions, reaction_names, lbs, ubs, gene_associations = sbml_content_handler.get_species_reactions(), sbml_content_handler.get_species_stoichiometries(), sbml_content_handler.get_num_species(), sbml_content_handler.get_num_reactions(), sbml_content_handler.get_reaction_names(), sbml_content_handler.get_lbs(), sbml_content_handler.get_ubs(), sbml_content_handler.get_gene_associations()

    return species_reactions, species_stoichiometries, num_species, num_reactions, reaction_names, lbs, ubs, gene_associations

def fetchSmatrix(species_reactions,species_stoichiometries,num_reactions,num_species):
	#Returns sparse S matrix from SBML in scipy dok format
	Sdok = dok_matrix((num_species,num_reactions), dtype=np.float32)
	for i in range(len(species_reactions)):
		row = species_reactions[i]
		k = 0
		for j in row:
			Sdok[i,j] = species_stoichiometries[i][k]
			k += 1

	#Convert DOK to CSR
	Scsr = Sdok.tocsr()

	return(Scsr)

def retrieve_information_from_SBML(model,default_bound=1000000):
	species_reactions, species_stoichiometries, num_species, num_reactions, reaction_names, lbs, ubs, gene_associations = _initialise_model(model,default_bound)
	Smatrix = fetchSmatrix(species_reactions,species_stoichiometries,num_reactions,num_species)
	
	return(Smatrix, num_species, num_reactions, reaction_names, lbs, ubs, gene_associations)
