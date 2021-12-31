# get_annotation.py

import collections
import re


ObjectAnnotation = collections.namedtuple('ObjectAnnotation',
                                        ['id', 'object_type', 'annotation'],
                                        )
def getSBO(sbo_num):
  """
  Reformat an SBO term into str.
  Return None if -1 (not provided).

  Parameters
  ----------
  sbo_num: int

  Returns
  -------
  list-tuple: [(object_id, 'sbo', SBO id)]
      Return [] if sbo term is -1 
      (i.e.,  not provided)
  """
  if sbo_num == -1:
    return []
  else:
    return [('sbo', 'SBO:' + format(sbo_num, '07d'))]
    
def getOntologyFromString(string_annotation):
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

def getCHEBIFromString(input_str):
  """
  Parses string and returns a CHEBI identifier. 
  If not, return None

  Parameters
  ----------
  str: string_annotation

  Returns
  -------
  str (CHEBI Id)
      Return None if none is provided
  """
  ontologies = getOntologyFromString(input_str)
  chebi_list = [val for val in ontologies if val[0]=='chebi']
  if chebi_list:
    return chebi_list[0][1]
  else:
    return None

def getQualifierFromString(input_str, qualifier):
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
  ontologies = getOntologyFromString(input_str)
  qualifier_list = [val for val in ontologies if val[0]==qualifier]
  if qualifier_list:
    return [val[1] for val in qualifier_list]
  else:
    return None