# test_reaction_recommender.py
# test for reaction_recommender.py 

import libsbml
import os
import unittest
from reaction_recommender import annotation_container as ac
from reaction_recommender import reaction_recommender as recom


BIGG_ECOLI = "e_coli_core.xml"

class TestReactionRecommender(unittest.TestCase):

  def setUp(self):
    reader = libsbml.SBMLReader()
    document = reader.readSBML(BIGG_ECOLI)
    self.model = document.getModel()
    # self.container = ac.AnnotationContainer(model_file=BIGG_ECOLI)
    self.raw_str = self.model.getSpecies('M_glc__D_e').getAnnotationString()
    self.recommender = recom.ReactionRecommender(model_file=BIGG_ECOLI)
    self.null_recommender = recom.ReactionRecommender()

  def testGetReferenceMatrix(self):
    dummy_dict = {'RHEA:10003': ['N', 'O', 'C5NO', 'C5O2'],
                  'RHEA:10007': ['C8NS', 'C8NS']}
    dummy_mat = self.recommender.getReferenceMatrix(dummy_dict)
    self.assertEqual(dummy_mat.shape, (2, 5))
    self.assertEqual(dummy_mat.loc['RHEA:10003', 'N'], 1)
    self.assertEqual(dummy_mat.loc['RHEA:10007', 'C8NS'], 1)

  def testGetCHEBIFormula(self):
  	one_formula = self.recommender.getCHEBIToFormula(id_to_chebis={'a': ['CHEBI:15379']})
  	self.assertEqual(one_formula['a'], ['O2'])
  	orig_formula = self.recommender.getCHEBIToFormula(id_to_chebis=None)
  	self.assertEqual(orig_formula['M_o2_c'], ['O2'])

  def testGetQueryMatrix(self):
    dummy_dict = {'RHEA:10003': ['N', 'O', 'C5NO', 'C5O2'],
                  'RHEA:10007': ['C8NS', 'C8NS']}
    dummy_mat = self.recommender.getReferenceMatrix(dummy_dict)
    dummy_query = self.recommender.getQueryMatrix(dummy_mat)
    self.assertEqual(dummy_query.shape, (3790, 5))

  def testGetCandidatesByReactionId(self):
    one_reac_cands = self.recommender.getCandidatesByReactionId('R_PFK')
    self.assertTrue(isinstance(one_reac_cands, dict))
    self.assertEqual(len(one_reac_cands['R_PFK']), 5)
    two_reac_cands = self.recommender.getCandidatesByReactionId(['R_PFK', 'R_PFL'])
    self.assertEqual(len(two_reac_cands), 2)
    all_reac_cands = self.recommender.getCandidatesByReactionId()
    self.assertEqual(len(all_reac_cands), 95)

  def testGetCandidateReport(self):
  	one_dict = {'R_PFK': ['RHEA:12423',
  	                      'RHEA:13380']}
  	one_rep = self.recommender.getCandidateReport(one_dict)
  	report = 'Reaction R_PFK: M_atp_c + M_f6p_c = M_adp_c + M_fdp_c + M_h_c\n'
  	report = report + 'Has possible RHEA IDs as below:\n\n0. '
  	report = report + '<RHEA:12420> ATP + D-tagatofuranose 6-phosphate = ADP + D-tagatofuranose 1,6-bisphosphate + H(+)\n'
  	report = report + '1. <RHEA:13377> alpha-D-glucose 1-phosphate + ATP = ADP + alpha-D-glucose 1,6-bisphosphate + H(+)\n'
  	report = report + '*********************************************************************************************************************\n'
  	self.assertEqual(one_rep, report)

  def testSortCandidates(self):
    sorted_cands = self.recommender.sortCandidates(['RHEA:10015', 'RHEA:10003'])
    self.assertEqual(sorted_cands, ['RHEA:10003', 'RHEA:10015'])

  def testFormatElement(self):
    one_reaction = self.model.getReaction('R_ATPS4r')
    adp = one_reaction.getReactant('M_adp_c')
    self.assertEqual(self.recommender.formatElement(adp),
                    'M_adp_c')
    hydro = one_reaction.getReactant('M_h_e')
    self.assertEqual(self.recommender.formatElement(hydro),
                    '4.0 * M_h_e')

  def testGetReactionString(self):
    one_reaction = self.model.getReaction('R_PGI')
    reaction_string = self.recommender.getReactionString(one_reaction)
    self.assertEqual(reaction_string, 'M_g6p_c = M_f6p_c')

if __name__ == '__main__':
  unittest.main()















