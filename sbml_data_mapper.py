'''
Requires:
path to SBML model, path to genes file (see format below)

Import SBML file and return all lists and matrices (including S matrix) required for constraint-based modelling
Import gene expression data from genes_file, a tab-delimited file of three columns with headings:
gene	mean	std

Merge these such that a new lists are created of size num_reactions (rxn_exp, rxn_exp_sds)

Returns:
Smatrix (sparse CSR), num_species (int), num_reactions (int), reaction_names (list), lbs (list), ubs (list), rxn_mean_expression (list), rxn_stdev_expression (list)
'''

import csv
import math

import SBMLRead as sbml_read
import gene_reaction_mapper


def _read_and_map_sbml_and_expression_data(model, genes_file, scaleName, scalingFactor, default_bound=1000000):
    '''
    model = sbml file
    genes_file = data file in format: gene	mean	std
    scaleName = name of reaction to scale to
    scalingFactor = multiplier for scaling reaction
    '''
    Smatrix, num_species, num_reactions, reaction_names, lbs, ubs, gene_associations = sbml_read.retrieve_information_from_SBML(
        model, default_bound)
    rxn_mean_expression, rxn_stdev_expression = _read_gene_expression_data(
        genes_file, gene_associations)

    # Rescale gene expression data:
    reaction_index_to_scale = reaction_names.index(scaleName)
    if math.isnan(rxn_mean_expression[reaction_index_to_scale]):
        print '\tWarning', reaction_index_to_scale, 'is NaN. All target vector components also will become NaN. DO NOT DO THIS!'
    rxn_mean_expression = [(x / rxn_mean_expression[reaction_index_to_scale])
                           * scalingFactor for x in rxn_mean_expression]
    rxn_stdev_expression = [(x / rxn_stdev_expression[reaction_index_to_scale])
                            for x in rxn_stdev_expression]  # This does nothing!

    return(Smatrix, num_species, num_reactions, reaction_names, lbs, ubs, rxn_mean_expression, rxn_stdev_expression)


def _read_and_map_sbml_and_relative_data(model, genes_file, default_bound=1000000):
    '''
    model = sbml file
    genes_file = data file in format: gene	mean	std
    '''
    Smatrix, num_species, num_reactions, reaction_names, lbs, ubs, gene_associations = sbml_read.retrieve_information_from_SBML(
        model, default_bound)
    rxn_ratio = _read_gene_relative_data(genes_file, gene_associations)

    return(Smatrix, num_species, num_reactions, reaction_names, lbs, ubs, rxn_ratio)


def _read_gene_expression_data(genes_file, gene_associations):

    gene_to_expression = {}
    gene_to_expression_sd = {}

    with open(genes_file, 'rU') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            gene_name = row['gene'].replace('-', '_')
            gene_to_expression[gene_name] = float(row['mean'])
            gene_to_expression_sd[gene_name] = float(row['std'])

        rxn_mean_expression, rxn_stdev_expression = gene_reaction_mapper.map_genes_to_reactions(
            gene_associations, gene_to_expression, gene_to_expression_sd)

    return(rxn_mean_expression, rxn_stdev_expression)


def _read_gene_relative_data(genes_file, gene_associations):

    gene_to_expression = {}
    gene_to_expression_sd = {}

    with open(genes_file, 'rU') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            gene_name = row['gene'].replace('-', '_')
            gene_to_expression[gene_name] = float(row['ratio'])
            gene_to_expression_sd[gene_name] = float(0)

        rxn_ratio, _ = gene_reaction_mapper.map_genes_to_reactions(
            gene_associations, gene_to_expression, gene_to_expression_sd)

    return(rxn_ratio)
