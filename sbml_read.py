'''
Fetches information from SBML by;

Smatrix, num_species, num_reactions, reaction_names, lbs, ubs,
gene_associations = retrieve_information_from_SBML_initialise_model(
path_to_sbml_file,default_bound)

e.g.
Smatrix, num_species, num_reactions, reaction_names, lbs, ubs,
gene_associations = retrieve_information_from_SBML_initialise_model(
'../Dave/mcisb-constraint-based-modelling/data/yeast/yeast_5.21_MCISB.xml',
100000000)

Creates:
  S matrix = a stoichiometric matrix in sparse CSR format of num_species rows
  and num_reactions columns
  num_species = the number of species in the SBML
  num_reactions = the number of reactions in the SBML
  reaction_names = a list of length num_reactions containing the names of
  each reaction in the SBML
  lbs, ubs = lists of length num_reactions containing the lower and upper
  bounds for each reaction
  gene_associations = a list of length num_reactions containing the genes
  associated with each reaction
'''
import xml.sax

from scipy.sparse import dok_matrix

import numpy as np
from sbml_content_handler import SBMLContentHandler


def retrieve_information_from_sbml(model, default_bound=10000):
    '''Parse SBML.'''
    species_reactions, species_stoichiometries, num_species, num_reactions, \
        reaction_names, lbs, ubs, gene_associations = _initialise_model(
            model, default_bound)

    s_matrix = _fetch_s_matrix(
        species_reactions, species_stoichiometries, num_reactions, num_species)

    return s_matrix, num_species, num_reactions, reaction_names, lbs, ubs, \
        gene_associations


def _initialise_model(sbml, default_bound=10000):
    handler = SBMLContentHandler(default_bound)
    xml.sax.parse(open(sbml), handler)

    return handler.get_species_reactions(), \
        handler.get_species_stoichiometries(), \
        handler.get_num_species(), \
        handler.get_num_reactions(), \
        handler.get_reaction_names(), \
        handler.get_lbs(), \
        handler.get_ubs(), \
        handler.get_gene_associations()


def _fetch_s_matrix(species_reactions, species_stoichiometries, num_reactions,
                    num_species):
    '''Returns sparse S matrix from SBML in scipy dok format.'''
    s_dok = dok_matrix((num_species, num_reactions), dtype=np.float32)

    for i, row in enumerate(species_reactions):
        k = 0

        for j in row:
            s_dok[i, j] = species_stoichiometries[i][k]
            k += 1

    # Convert DOK to CSR
    return s_dok.tocsr()
