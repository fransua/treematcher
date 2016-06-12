from __future__ import absolute_import
from ete3 import Tree, PhyloTree, NCBITaxa
from treematcher import TreePattern
import unittest
import sys

# DATABASE_PATH = "testdb.sqlite"
DATABASE_PATH = "../ete/testdb.sqlite"

class Test_TreePattern(unittest.TestCase):
    def test_ete_params(self):
        """
        tests exact match is working
        tests basic ete parameters like node.name, node.dist, node.children are working

        """
        pattern0 = """
            exact ( bye , kk );
            """
        pattern1 = """
            exact ( hello , kk );
            """
        pattern2 = """
            (
            'len(@.children) > 2 and @.name in ("hello","kk") '
            )
            '(len(@.name) > 3) and @.dist >= 0.5';
            """

        pattern0 = TreePattern(pattern0, format=8, quoted_node_names=False)
        pattern1 = TreePattern(pattern1, format=8, quoted_node_names=False)
        pattern2 = TreePattern(pattern2, format=8, quoted_node_names=True)

        tree = Tree("(hello,(1,2,3)kk)pasa:1;", format=1)
        self.assertEqual(pattern0.find_match(tree, None)[0], False)
        self.assertEqual(pattern1.find_match(tree, None)[0], True)
        self.assertEqual(pattern2.find_match(tree, None)[0], True)

    def test_syntax_to_python(self):
        """
        tests syntax we've created like node.leaves, node.children are working
         """

        pattern1 = """( '  @.dist >= 0.5 ' , ' @.dist<2  ')
            '    "Pan_troglodytes_1" in @.leaves and "Homo_sapiens_1" in @.children ';"""

        pattern2 = """( '  Distance greater than or equal to 0.5 ' , ' Distance less than 2 ' )
            '     "Pan_troglodytes_1" in Leaves and "Homo_sapiens_1" in Children';"""

        pattern1 = TreePattern(pattern1, format=8, quoted_node_names=True)
        pattern2 = TreePattern(pattern2, format=8, quoted_node_names=True)

        tree = PhyloTree(
            "((((Anolis_carolinensis_1:1, Gallus_gallus_1:1), (Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_1:1))), ((Danio_rerio_1:1, (Xenopus_laevis_1:1, Anolis_carolinensis_1:1)), Saccharomyces_cerevisiae_2:1)), Saccharomyces_cerevisiae_1:1);",
            format=1)
        self.assertEqual(pattern1.find_match(tree, None)[0], True)
        self.assertEqual(pattern2.find_match(tree, None)[0], True)
        tree = PhyloTree(
            "((((Anolis_carolinensis_1:1, Gallus_gallus_1:1), (Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_2:1))), ((Danio_rerio_1:1, (Xenopus_laevis_1:1, Anolis_carolinensis_1:1)), Saccharomyces_cerevisiae_2:1)), Saccharomyces_cerevisiae_1:1);",
            format=1)
        self.assertEqual(pattern1.find_match(tree, None)[0], False)
        self.assertEqual(pattern2.find_match(tree, None)[0], False)

    def test_species_taxonomy(self):
        """
        tests if node.species and ncbi_query are working
        """

        # test node.species

        species_tree = PhyloTree(
            "((Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_1:1), Saccharomyces_cerevisiae_1:1));",
            format=1)
        species_tree.set_species_naming_function(lambda n: n.name.split("_")[1] if "_" in n.name else '')

        pattern0 = """((( ' @.species in ("sapiens","pygmaeus")  '))' "Pan_troglodytes_1" ');"""
        pattern0 = TreePattern(pattern0, format=8, quoted_node_names=True)

        self.assertEqual(pattern0.find_match(species_tree, None)[0], True)

        # test ncbi taxonomy

        ncbi = NCBITaxa(dbfile=DATABASE_PATH)
        taxonomy_tree = PhyloTree("((9598, 9606), 10090);", sp_naming_function=lambda name: name)
        taxonomy_tree.annotate_ncbi_taxa(dbfile=DATABASE_PATH)
        root = taxonomy_tree.get_tree_root()

        pattern1 = """ '  @.sci_name == "Euarchontoglires" ';"""
        pattern2 = """
          (( ' @.sci_name=="Homo sapiens" , 9526 in @.lineage ' )' @.rank=="subfamily" and @.taxid == 207598 ')
          '  @.sci_name == "Euarchontoglires" and "cellular organisms" in @.named_lineage';
          """

        pattern1 = TreePattern(pattern1, format=1, quoted_node_names=True)
        pattern2 = TreePattern(pattern2, format=1, quoted_node_names=True)
        match1 = pattern1.find_match(taxonomy_tree, None)
        match2 = pattern2.find_match(taxonomy_tree, None)

        self.assertEqual(match1, (True, root))
        self.assertEqual(match2, (True, root))

    def test_custom_fuctions(self):
        """
            tests custom functions are working
        """
        custom_functions = {"length": length}

        pattern = """
            (
            'len(@.children) > 2 and @.name in ("hello","bye") '
            )
            '(length(@.name) < 3) and @.dist >= 0.5';
            """

        pattern = TreePattern(pattern, format=8, quoted_node_names=True)
        tree = Tree("(hello,(1,2,3)kk)pasa:1;", format=1)
        self.assertEqual(pattern.find_match(tree, custom_functions)[0], False)

        tree = Tree("((kk,(1,2,3)bye)y:1, NODE);", format=1)
        self.assertEqual(pattern.find_match(tree, custom_functions)[0], True)

        tree = Tree("(((1,2,3)bye)y:1, NODE);", format=1)
        self.assertEqual(pattern.find_match(tree, custom_functions)[0], True)

        tree = Tree("(((1,2,3)bye,kk)y:1, NODE);", format=1)
        self.assertEqual(pattern.find_match(tree, custom_functions)[0], True)

    def test_evol_events(self):
        """
      tests that .evol_event returns "R" or "D"
      """
        pass

    def test_thousands_of_trees(self):

        """
        Search over 26 thousand trees (naming format: NumericTaxid.SequenceName)
        for nodes containing branches that separate two groups of primate genes where,
        in one side, the human gene has been lost,
        and the branch support value of the matching node is higher than 0.9.

                                  /-Any primate taxid (9443 in lineage)
        support >= 0.9--|
                                  \-Any primate taxid except human

        """
        try:
            sample_list = "1000Trees.tsv"  # raw data list from Alan's lab [id, family, genus, species]

        except IndexError:
            print "invalid filename"



        with open(sample_list, 'r') as sample:
            lines = sample.readlines()
            line_count = 0
            match_count = 0

            for line in lines:
                line_count = line_count + 1
                cells = line.split('\t')
                # project = cells[0]
                # id = cells[1]
                # software = cells[2]
                # trim = cells[3]
                tree = cells[4]

                tree = PhyloTree(tree, format=1, quoted_node_names=False)

                tree.set_species_naming_function(lambda n: n.name.split(".")[0] if "." in n.name else '')
                tree.annotate_ncbi_taxa(dbfile=DATABASE_PATH)

                pattern = """
                    ( ' 9443 in @.lineage ' , ' 9443 in @.lineage and @.name!=9606 ' )' @.support >= 0.9 ';
                    """
                pattern = TreePattern(pattern, format=8, quoted_node_names=False)

                if pattern.find_match(tree, None)[0]:
                    match_count = match_count + 1

                if (line_count % 200 == 0):
                    print str(match_count) + " matches found out of " + str(line_count)

            print "finished with " + str(match_count) + " matches for " + str(line_count) + " trees."

        self.assertEqual(match_count, line_count)



def length(txt):
    return len(txt)


def run():
    unittest.main()

if __name__ == '__main__':
    run()
