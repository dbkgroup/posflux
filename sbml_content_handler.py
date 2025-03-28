'''
(c) University of Liverpool 2019

All rights reserved.
'''
# pylint: disable=too-many-instance-attributes
import xml.sax


class SBMLContentHandler(xml.sax.ContentHandler):
    '''SBML SAX parser'''

    def __init__(self, _default_bound):
        xml.sax.ContentHandler.__init__(self)
        self.__default_bound = _default_bound
        self.__species_reactions = None
        self.__species_stoichiometries = None
        self.__num_species = 0
        self.__num_reactions = 0
        self.__reaction_names = []
        self.__lbs = []
        self.__ubs = []
        self.__gene_associations = []

        self.__species_ids = []
        self.__in_p = None
        self.__in_gene_association = False
        self.__in_reactants = None

        self.__gene_association = ''

    def startElement(self, name, attrs):
        if name == 'species' \
            and ('boundaryCondition' not in attrs
                 or not attrs.getValue('boundaryCondition')):
            self.__species_ids.append(attrs.getValue('id'))
        elif name == 'reaction':
            self.__reaction_names.append(attrs.getValue('name')
                                         if 'name' in attrs else '')
        elif name == 'p':
            self.__in_p = True
        elif name == 'listOfReactants':
            self.__in_reactants = True
        elif name == 'listOfProducts':
            self.__in_reactants = False
        elif name == 'speciesReference':
            species_id = attrs.getValue('species')

            if species_id in self.__species_ids:
                i = self.__species_ids.index(species_id)

                stoichiometry = float(attrs.getValue('stoichiometry')) \
                    if 'stoichiometry' in attrs else 1

                self.__species_reactions[i].append(self.__num_reactions)
                self.__species_stoichiometries[i].append(
                    -stoichiometry if self.__in_reactants else stoichiometry)

        elif name == 'parameter' and attrs.getValue('id') == 'LOWER_BOUND':
            self.__lbs.append(
                max(float(attrs.getValue('value')), -self.__default_bound))
        elif name == 'parameter' and attrs.getValue('id') == 'UPPER_BOUND':
            self.__ubs.append(
                min(float(attrs.getValue('value')), self.__default_bound))

    def characters(self, content):
        if self.__in_p and (content.startswith('GENE ASSOCIATION:')
                            or content.startswith('GENE_ASSOCIATION:')):
            self.__in_gene_association = True
        if self.__in_gene_association:
            self.__gene_association += content[len('GENE_ASSOCIATION:'):] \
                if 'GENE_ASSOCIATION:' in content \
                else content[len('GENE ASSOCIATION:'):] \
                if 'GENE ASSOCIATION:' in content \
                else content
            self.__gene_association = self.__gene_association.replace('-', '_')

    def endElement(self, name):
        if name == 'listOfSpecies':
            self.__num_species = len(self.__species_ids)
            self.__species_reactions = [[] for _ in range(self.__num_species)]
            self.__species_stoichiometries \
                = [[] for _ in range(self.__num_species)]
        elif name == 'p':
            self.__in_p = False

            if self.__in_gene_association:
                self.__in_gene_association = False
                self.__gene_associations.append(self.__gene_association)
                self.__gene_association = ''

        elif name == 'reaction':
            self.__num_reactions = self.__num_reactions + 1

    def get_species_reactions(self):
        '''TODO'''
        return self.__species_reactions

    def get_species_stoichiometries(self):
        '''TODO'''
        return self.__species_stoichiometries

    def get_num_species(self):
        '''TODO'''
        return self.__num_species

    def get_num_reactions(self):
        '''TODO'''
        return self.__num_reactions

    def get_reaction_names(self):
        '''TODO'''
        return self.__reaction_names

    def get_lbs(self):
        '''TODO'''
        return self.__lbs

    def get_ubs(self):
        '''TODO'''
        return self.__ubs

    def get_gene_associations(self):
        '''TODO'''
        return self.__gene_associations
