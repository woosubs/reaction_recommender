# reaction_recommender.py
# clsss for predicting reaction annotations

import itertools
import libsbml
import numpy as np
import os
import pandas as pd
import pickle

# below should exist as a separate class
import annotation_container as ac
import constants as cn
import get_annotation as ga


# Nov 14, 2021: using a shortened version, meaning H and D are removed
with open(os.path.join(cn.DATA_DIR, 'all_shortened_chebi_to_formula_19nov2021.pkl'), 'rb') as f:
  ref_shortened_chebi_to_formula = pickle.load(f)
# map secondary id to primary id
with open(os.path.join(cn.DATA_DIR, 'chebi_second2prime_8nov2021.pickle'), 'rb') as f:
  ref_second2prime_dict = pickle.load(f)
with open(os.path.join(cn.RHEA_DIR, 'rhea2chebi_reference.pkl'), 'rb') as f:
  ref_rhea2chebi = pickle.load(f)
with open(os.path.join(cn.RHEA_DIR,'rhea_all2master.pkl'), 'rb') as f:
  ref_rhea2master = pickle.load(f)

# below needs to be incorported with class
def getRheaToFormula(ref_rhea_to_chebi, ref_chebi_to_formula):
  """
  Get a dictionary of RHEA-components (formulae)
  using RHEA->CHEBI, and CHEBI->Formula map.
  :param dict ref_rhea_to_chebi:
  :param dict ref_chebi_to_formula:
  :return dict
  """
  rhea_to_formula = dict()
  for one_rhea in ref_rhea_to_chebi.keys():
    one_chebis = ref_rhea_to_chebi[one_rhea]
    one_formulas = [ref_chebi_to_formula[k] for k in one_chebis\
                    if k in ref_chebi_to_formula.keys() and ref_chebi_to_formula[k]]
    # remove None
    filt_one_formulas = [val for val in one_formulas if val]
    rhea_to_formula[one_rhea] = filt_one_formulas
  return rhea_to_formula

ref_short_rhea2formula = getRheaToFormula(ref_rhea_to_chebi=ref_rhea2chebi,
                                          ref_chebi_to_formula=ref_shortened_chebi_to_formula)


class ReactionRecommender(object):
  """
  Recommend reaction annotations
  using given list of species names.
  Can be used in two ways:
    1. use data from SBML model file
    2. user's input (list of species annotations)

  Attributes
  ----------
  species: list of species for annotation
  reactions: list of reactions for annotation
  spec_dict: dictionary, key=species ID, itm=CHEBI annotation
  spec_formula_dict: dictionary, key=species ID, itm=formula
  ref_mat: Reference matrix
  """

  def __init__(self,
               model_file=None):
    """
    Parameters
    ----------
    model_file: str
        Address/name of the .xml model file
    """
    # model file can take None
    if model_file is None:
      self.spec_dict = None
      self.reac_dict = None
      self.spec_formula_dict = None
    else:
      self.container = ac.AnnotationContainer(model_file)
      self.spec_dict = self.container.spec_dict
      self.reac_dict = self.container.reac_dict
      self.spec_formula_dict = self.getCHEBIToFormula()
    self.ref_mat = self.getReferenceMatrix()


    #   # TODO: some of below operation will be done
    #   # by the annotation container class
    #   reader = libsbml.SBMLReader()
    #   document = reader.readSBML(model_file)
    #   model = document.getModel()
    #   self.model = model
    #   spec_dict = dict()
    #   for one_spec in model.getListOfSpecies():
    #     spec_dict[one_spec.getId()] = ga.getQualifierFromString(one_spec.getAnnotationString(), 'chebi')
    #   self.spec_dict = spec_dict
    #   ###### below needed! ######
    #   self.spec_formula_dict = self.getCHEBIToFormula()
    #   ###### above needed! ### ###            
    #   reac_dict = dict()
    #   for one_reaction in model.getListOfReactions():
    #     reactants = [val.species for val in one_reaction.getListOfReactants()]
    #     products = [val.species for val in one_reaction.getListOfProducts()]
    #     reac_dict[one_reaction.getId()] = list(set(reactants + products))
    #   self.reac_dict = reac_dict
    # self.ref_mat = self.getReferenceMatrix()

  def getReferenceMatrix(self,
                         ref_reaction_to_components=ref_short_rhea2formula):
    """
    Create a matrix (DataFrame)
    with the relationship between
    reaction identifier to components.

    Parameters
    ----------
    ref_reaction_to_components: dict
        Dictionary that maps reaction identifiers
        to components (identififer or formula)

    Returns
    -------
    ref_df: pandas.DataFrame
        Matrix with ones/zeros to indicate
        the relationships
    """  	
    all_components = set()
    for one_id in ref_reaction_to_components.keys():
      all_components = all_components.union(set(ref_reaction_to_components[one_id]))

    ref_df = pd.DataFrame(0, 
                         index=ref_reaction_to_components.keys(),
                         columns=all_components)
    # fill the cells using ref_reaction_to_components  
    for one_id in ref_reaction_to_components.keys():
      one_components = ref_reaction_to_components[one_id]
      for val in one_components:
        if val in ref_df.columns:
          ref_df.loc[one_id, val] = 1
    return ref_df

  def getCHEBIToFormula(self, 
  	                    ref_chebi_secondary_to_primary=ref_second2prime_dict,
  	                    ref_chebi_to_formula=ref_shortened_chebi_to_formula,
  	                    id_to_chebis=None):
    """
    Create a dictionary, 
    mapping species id to corresponding formula.
    Relationships need to be defined by ref_ dictionaries.
    A previous relationship defining 
    Species ID (str, dictionary key) 
    to CHEBIs (str-list) should be given.

    Parameters
    ----------
    ref_chebi_secondary_to_primary: dict
        Dictionary of CHEBI terms,
        mapping all chebi terms to 
        their corresponding primary terms.

    ref_chebi_to_formula: dict
        Dictionary mapping primary chebi terms to
        their corresponding formula.

    id_to_chebis: dict
        Dictionary, mapping species id to CHEBI terms

    Returns
    -------
    id_to_formula: dict
        Dictionary mapping given species ID - formula.
    """
    if id_to_chebis is None:
      id_to_chebis = self.spec_dict
    filt_id_to_chebis = {one_k:id_to_chebis[one_k] for one_k in id_to_chebis.keys()\
                         if id_to_chebis[one_k]}
    none_id_to_chebis = {one_k:None for one_k in id_to_chebis.keys()\
                         if id_to_chebis[one_k] is None}
    # collect all chebi terms and create a dictionary
    all_chebis = set()
    for one_k in filt_id_to_chebis.keys():
      all_chebis = all_chebis.union(set(filt_id_to_chebis[one_k]))
    all_chebi2formula = {val:ref_chebi_to_formula[ref_chebi_secondary_to_primary[val]] for val in all_chebis\
                         if val in ref_chebi_secondary_to_primary.keys()}

    # dictionary: species ID to [chebi -> formula mapped] list
    id_to_formula = {val:list({all_chebi2formula[one_chebi] for one_chebi in filt_id_to_chebis[val]\
    	             if one_chebi in all_chebi2formula.keys()})\
                     for val in filt_id_to_chebis.keys()}
    id_to_formula.update(none_id_to_chebis)
    return id_to_formula

  def getQueryMatrix(self,
  	                 reac_comp_dict,
  	                 ref_mat=None):
    """
    Get a query matrix for calculation. 
    If species exists for a given reaction,
    corresponding cell has 1; otherwise 0.

    Parameters
    ----------
    reac_comp_dict: dict
        Dictionary of reaction (key) and
        species/formula components (item).
        Can be obtained from getReactionToSpeciesRelationship().

    Returns
    -------
    query_mat: pandas.DataFrame/None
               Query matrix to multiply to the reference matrix
    """
    if ref_mat is None:
      ref_mat = self.ref_mat
    if ref_mat is not None:
      query_mat = pd.DataFrame(0,
                               index=ref_mat.columns,
                               columns=reac_comp_dict.keys())
      for k in reac_comp_dict.keys():
        items = [val for val in reac_comp_dict[k] if val in query_mat.index]
        for one_itm in items:
          query_mat.loc[one_itm, k] = 1
    else:
      query_mat = None
    return query_mat

  def getCandidatesByReactionId(self,
  	                            reaction_ids=None):
    """
    Get the candidates of given reaction IDs.
    IDs are from the SBML model

    Parameters
    ----------
    reactions: str/str-list
        List of reaction IDs of a model.
        If nothing is provided, 
        it will use all reactions.

    Returns
    -------
    candidates: dict
        Ordered list of candidates..
    """
    if reaction_ids is None:
      reaction_ids = list(self.reac_dict.keys())
    elif isinstance(reaction_ids, str):
      reaction_ids = [reaction_ids]
    sub_reac_dict = {val:self.reac_dict[val] for val in reaction_ids}
    query_reac_to_formula = {val:list(itertools.chain(*[self.spec_formula_dict[itm] for itm in sub_reac_dict[val] if self.spec_formula_dict[itm]]))\
                             for val in sub_reac_dict.keys()}
    query_df = self.getQueryMatrix(reac_comp_dict=query_reac_to_formula,
    	                           ref_mat=self.ref_mat)
    # multiply two matrices
    multi_mat = self.ref_mat.dot(query_df)
    maxes = multi_mat.max()
    candidates_dict = dict()
    for one_idx in maxes.index:
      one_multi = multi_mat.loc[:,one_idx]
      candidates = one_multi[one_multi==maxes[one_idx]].index
      candidates_dict[one_idx] = list(candidates)
    return candidates_dict

  def getCandidateReport(self, candidates,
  	                     ref_rhea_df,
  	                     ref_rhea_to_master=ref_rhea2master):
    """
    Get reports for candidates.

    Parameters
    ----------
    candidates: dict
        Dictionary mapping reaction IDs
        to a list of RHEA identifiers.
    ref_rhea_to_master: dict
        Dictionary mapping RHEA IDs
        to its corresponding master(default) ID.

    Returns
    -------
    report: str
        Report of candidates
    """
    report = ""
    for idx, one_k in enumerate(candidates.keys()):
      if idx > 0:
        report = report + "\n"
      report = report + "Reaction %s: " % one_k
      report = report + self.getReactionString(self.container.model.getReaction(one_k))
      sorted_candidates = self.sortCandidates(candidate_list=candidates[one_k])
      one_rheas = [ref_rhea_to_master[val] for val in sorted_candidates]
      report = report + "\nHas possible RHEA IDs as below:\n\n"
      for idx, one_rhea in enumerate(one_rheas):
      	if one_rhea in ref_rhea_df.index:
          report = report + "%d. <%s> %s\n" % (idx, one_rhea, ref_rhea_df.loc[one_rhea, 'equation'])
      report = report + "***************************************"*3 + "\n"
    return report

  def sortCandidates(self, candidate_list):
    """
    Sort candidate list (RHEA IDs)
    using ascending 
    1. number of nonzero cells
    2. RHEA ID number

    Parameters
    ----------
    candidate_list: str-list 
                    e.g., ['RHEA:12345', 'RHEA:23456']
    """
    num_ele = [len(np.nonzero(np.array(self.ref_mat.loc[val,:]))[0]) for val in candidate_list]
    num_rhea = [int(val[5:]) for val in candidate_list]
    return [candidate_list[idx] for idx in np.lexsort((num_rhea, num_ele))]

  def formatElement(self, one_species):
    """
    Formst species stoichiomety and species
    """
    element_stoichiometry = one_species.stoichiometry
    if element_stoichiometry==1.0:
      return one_species.species
    else:
      return "%s * %s" % (element_stoichiometry, one_species.species)

  def getReactionString(self, one_reaction):
    """
    Get reaction as a string
    to print out
    """
    reactants = " + ".join([self.formatElement(spec) for spec in one_reaction.getListOfReactants()])
    products = " + ".join([self.formatElement(spec) for spec in one_reaction.getListOfProducts()])
    result = reactants + " = " + products
    return result

