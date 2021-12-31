# annotation_container.py
# Container to parse and store SBML annotations

import collections
import libsbml
import re


class AnnotationContainer(object):
  """
  Parses an SBML model (in XML),
  and stores annotations

  Attributes
  ----------
  species: list of species for annotation (in CHEBI)
  reactions: list of reactions for annotation
  spec_dict: dictionary, key=species ID, itm=CHEBI annotation
  """

  def __init__(self,
               model_file=None):
    """
    Parameters
    ----------
    model_file: str
        Address/name of the .xml model file
    """
    reader = libsbml.SBMLReader()
    document = reader.readSBML(model_file)
    model = document.getModel()
    self.model = model
    spec_dict = dict()
    for one_spec in model.getListOfSpecies():
      spec_dict[one_spec.getId()] = self.getQualifierFromString(one_spec.getAnnotationString(), 'chebi')
    self.spec_dict = spec_dict
    reac_dict = dict()
    for one_reaction in model.getListOfReactions():
      reactants = [val.species for val in one_reaction.getListOfReactants()]
      products = [val.species for val in one_reaction.getListOfProducts()]
      reac_dict[one_reaction.getId()] = list(set(reactants + products))
    self.reac_dict = reac_dict

  def getQualifierFromString(self, input_str, qualifier):
    """
    Parses string and returns a CHEBI identifier. 
    If not, return None

    Parameters
    ----------
    str: string_annotation

    Returns
    -------
    str (ontology Id)
        Return None if none is provided
    """
    ontologies = self.getOntologyFromString(input_str)
    qualifier_list = [val for val in ontologies if val[0]==qualifier]
    if qualifier_list:
      return [val[1] for val in qualifier_list]
    else:
      return None

  def getOntologyFromString(self, string_annotation):

    """
    Parse string and return string annotation,
    marked as <bqbiol:is> or <bqbiol:isVersionOf>.
    If neither exists, return None.

    Parameters
    ----------
    str: string_annotation

    Returns
    -------
    list-tuple (ontology type, ontology id)
         Return [] if none is provided
    """
    # first, extracts strings tagged as bqbiol:is or bqbiol:isVersionOf.
    is_str = ''
    isVersionOf_str = ''
    is_str_match = re.findall('<bqbiol:is[^a-zA-Z].*?<\/bqbiol:is>',
                              string_annotation,
                              flags=re.DOTALL)
    if len(is_str_match)>0:
      is_str_match_filt = [s.replace("      ", "") for s in is_str_match]
      is_str = '\n'.join(is_str_match_filt)

    is_VersionOf_str_match = re.findall('<bqbiol:isVersionOf[^a-zA-Z].*?<\/bqbiol:isVersionOf>',
                                        string_annotation,
                                        flags=re.DOTALL)
    #
    if len(is_VersionOf_str_match) > 0:
      is_VersionOf_str_match_filt = [s.replace("      ", "") for s in is_VersionOf_str_match]
      isVersionOf_str = '\n'.join(is_VersionOf_str_match_filt)
    #
    combined_str = is_str + isVersionOf_str
    if combined_str == '':
      return []
    identifiers_list = re.findall('identifiers\.org/.*/', combined_str)
    return [(r.split('/')[1],r.split('/')[2].replace('\"', '')) \
            for r in identifiers_list]


