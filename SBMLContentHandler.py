'''
Created on 23 Jun 2014

@author: neilswainston
'''
import xml.sax


class SBMLContentHandler(xml.sax.ContentHandler):
    '''SBML SAX parser'''
    _species_reactions = None
    _species_stoichiometries = None
    _num_reactions = 0
    _reaction_names = []
    _lbs = []
    _ubs = []
    _gene_associations = []

    _default_bound = None
    _species_ids = []
    _in_p = None
    _in_gene_association = False
    _in_reactants = None

    _gene_association = ''

    def __init__(self, _default_bound):
        xml.sax.ContentHandler.__init__(self)
        self._default_bound = _default_bound

    def startElement(self, name, attrs):
        if name == 'species' \
            and (not 'boundaryCondition' in attrs
                 or attrs.getValue('boundaryCondition') == 'false'):
            self._species_ids.append(attrs.getValue('id'))
        elif name == 'reaction':
            self._reaction_names.append(attrs.getValue('name')
                                        if 'name' in attrs else '')
        elif name == 'p':
            self._in_p = True
        elif name == 'listOfReactants':
            self._in_reactants = True
        elif name == 'listOfProducts':
            self._in_reactants = False
        elif name == 'speciesReference':
            species_id = attrs.getValue('species')

            if species_id in self._species_ids:
                i = self._species_ids.index(species_id)

                stoichiometry = float(attrs.getValue('stoichiometry')) \
                    if 'stoichiometry' in attrs else 1

                self._species_reactions[i].append(self._num_reactions)
                self._species_stoichiometries[i].append(
                    -stoichiometry if self._in_reactants else stoichiometry)

        elif name == 'parameter' and attrs.getValue('id') == 'LOWER_BOUND':
            self._lbs.append(
                max(float(attrs.getValue('value')), -self._default_bound))
        elif name == 'parameter' and attrs.getValue('id') == 'UPPER_BOUND':
            self._ubs.append(
                min(float(attrs.getValue('value')), self._default_bound))

    def characters(self, content):
        if self._in_p and (content.startswith('GENE ASSOCIATION:') or content.startswith('GENE_ASSOCIATION:')):
            self._in_gene_association = True
        if self._in_gene_association:
            self._gene_association += content[len('GENE_ASSOCIATION:'):] if 'GENE_ASSOCIATION:' in content \
                                else content[len('GENE ASSOCIATION:'):] if 'GENE ASSOCIATION:' in content \
                                else content
            self._gene_association = self._gene_association.replace('-', '_')

    def endElement(self, name):
        if name == 'listOfSpecies':
            self._num_species = len(self._species_ids)
            self._species_reactions = [[] for _ in range(self._num_species)]
            self._species_stoichiometries \
                = [[] for _ in range(self._num_species)]
        elif name == 'p':
            self._in_p = False

            if self._in_gene_association:
                self._in_gene_association = False
                self._gene_associations.append(self._gene_association)
                self._gene_association = ''

        elif name == 'reaction':
            self._num_reactions = self._num_reactions + 1

    def get_species_reactions(self):
        '''TODO'''
        return self._species_reactions

    def get_species_stoichiometries(self):
        '''TODO'''
        return self._species_stoichiometries

    def get_num_species(self):
        '''TODO'''
        return self._num_species

    def get_num_reactions(self):
        '''TODO'''
        return self._num_reactions

    def get_reaction_names(self):
        '''TODO'''
        return self._reaction_names

    def get_lbs(self):
        '''TODO'''
        return self._lbs

    def get_ubs(self):
        '''TODO'''
        return self._ubs

    def get_gene_associations(self):
        '''TODO'''
        return self._gene_associations
