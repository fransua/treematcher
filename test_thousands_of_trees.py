from __future__ import absolute_import
from ete3 import Tree, PhyloTree, NCBITaxa
from treematcher import TreePattern
import unittest
import sys

# DATABASE_PATH = "testdb.sqlite"
DATABASE_PATH = "../ete/testdb.sqlite"

class Test_TreePattern(unittest.TestCase):

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
            sample_list = "maNOG.trees.tsv"  # raw data list from Alan's lab [id, family, genus, species]

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
                del tree

            print "finished with " + str(match_count) + " matches for " + str(line_count) + " trees."

        self.assertEqual(match_count, line_count)


def length(txt):
    return len(txt)


def run():
    unittest.main()

if __name__ == '__main__':
    run()
