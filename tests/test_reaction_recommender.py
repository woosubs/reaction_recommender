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

if __name__ == '__main__':
  unittest.main()