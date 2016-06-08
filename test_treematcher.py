from __future__ import absolute_import
from ete3 import PhyloTree, NCBITaxa
from ete3.ncbi_taxonomy import ncbiquery
from treematcher import TreePattern
import unittest


DATABASE_PATH = "testdb.sqlite"


class Test_TreePattern(unittest.TestCase):

  def test_ete_params(self):
      """
      tests basic ete parameters like node.name, node.dist, node.children are working
      """
      pass

  def test_syntax_to_python(self):
      """
      tests syntax we've created like node.leaves, node.children are working
      """
      pass

  def test_species_taxonomy(self):
      """
      tests if node.species and ncbi_query are working
      """
      species_tree = PhyloTree("((((Anolis_carolinensis_1:1, Gallus_gallus_1:1), (Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_1:1))), ((Danio_rerio_1:1, (Xenopus_laevis_1:1, Anolis_carolinensis_1:1)), Saccharomyces_cerevisiae_2:1)), Saccharomyces_cerevisiae_1:1);", format=1)
      species_tree.set_species_naming_function(lambda n: n.name.split("_")[1] if "_" in n.name else '')


      pattern0 = """
          (
          (
          ( ' @.species in ("sapiens","pygmaeus")  ')
          )
          ' "Pan_troglodytes_1" '
          )
          ;
          """
      pattern0 = TreePattern(pattern0, format=8, quoted_node_names=True)

      self.assertEqual(pattern0.find_match(species_tree, None)[0], True)

      ncbi = NCBITaxa(dbfile=DATABASE_PATH)

      taxonomy_tree = PhyloTree("((9598, 9606), 10090);", sp_naming_function=lambda name: name)
      taxonomy_tree.annotate_ncbi_taxa(dbfile=DATABASE_PATH)
      root = taxonomy_tree.get_tree_root()


      pattern1 = """
          '  @.sci_name == "Euarchontoglires" '
          ;
          """

      pattern2 = """
          (
          ( ' @.sci_name=="Homo sapiens" , 9526 in @.lineage ' )
          ' @.rank=="subfamily" and @.taxid == 207598 '
          )
          '  @.sci_name == "Euarchontoglires" and "cellular organisms" in @.named_lineage'
          ;
          """

      pattern1 = TreePattern(pattern1, format=1, quoted_node_names=True)
      pattern2 = TreePattern(pattern2, format=1, quoted_node_names=True)
      match1 = pattern1.find_match(taxonomy_tree, None)
      match2 = pattern2.find_match(taxonomy_tree, None)
      self.assertEqual(match1, (True, root))
      self.assertEqual(match2, (True, root))

  def test_custom_fuctions(self):
      """
      tests if built-in custom functions are working
      """
      pass

  def test_evol_events(self):
      """
      tests that .evol_event returns "R" or "D"
      """
      pass


def run():
    unittest.main()

if __name__ == '__main__':
    run()
