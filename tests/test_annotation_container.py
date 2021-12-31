# test_annotation_container.py
# testing annotation_container.py

import libsbml
import os
import unittest
from reaction_recommender import annotation_container as ac

BIGG_ECOLI = "e_coli_core.xml"

class TestAnnotationContainer(unittest.TestCase):

  def setUp(self):
    reader = libsbml.SBMLReader()
    document = reader.readSBML(BIGG_ECOLI)
    self.model = document.getModel()
    self.container = ac.AnnotationContainer(model_file=BIGG_ECOLI)
    self.raw_str = self.model.getSpecies('M_glc__D_e').getAnnotationString()

  def testGetOntologyFromString(self):
    parsed_str = self.container.getOntologyFromString(self.raw_str)
    self.assertTrue(parsed_str)
    bigg_metas = [val for val in parsed_str if val[0]=='bigg.metabolite']
    self.assertEqual(bigg_metas[0][1], 'glc__D')
    self.assertEqual(len(parsed_str), 16)
    chebis = set([val[1] for val in parsed_str if val[0]=='chebi'])
    self.assertEqual(chebis, set(['CHEBI:12965', 'CHEBI:20999',
                                  'CHEBI:4167', 'CHEBI:17634'],
                                 ))

  def testGetQualifierFromString(self):
    bigg_metas_ele = self.container.getQualifierFromString(self.raw_str, 'bigg.metabolite')
    self.assertEqual(bigg_metas_ele, ['glc__D'])
    chebi_ele = self.container.getQualifierFromString(self.raw_str, 'chebi')
    self.assertEqual(set(chebi_ele), set(['CHEBI:12965', 'CHEBI:20999',
                                          'CHEBI:4167', 'CHEBI:17634'],
                                         ))


if __name__ == '__main__':
  unittest.main()
