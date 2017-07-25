import unittest

#assuming ete3 is not installed. they will be imported from local directory on developement modules.
from ete3 import PhyloTree, Tree, NCBITaxa
from treematcher import TreePattern

class Relax_number_of_children(unittest.TestCase):
    def setUp(self):
        tree = Tree(""" ((((A, B, C, D, E, F, G, H, I, J)), (I, J, K, L, M))) ; """)

        (tree&'A').dist = 0.5
        (tree&'B').dist = 0.5
        (tree&'C').dist = 0.5
        (tree&'D').dist = 0.5
        (tree&'E').dist = 0.3
        (tree&'F').dist = 0.3
        (tree&'G').dist = 0.1
        (tree&'H').dist = 0.1

        self.tree = tree

        self.pt1 = TreePattern(""" ( '@.name == "A", @.dist==0.5'  ); """, quoted_node_names=True)
        self.pt2 = TreePattern(""" ( '@.name == "A", @.dist==1'  ); """, quoted_node_names=True)
        self.pt3 = TreePattern(""" ( '@.name == "A", @.dist==0.5', '@.dist==1.1, +'); """, quoted_node_names=True)
        self.pt4 = TreePattern(""" ( '@.name == "A", @.dist==0.5', '@.dist==0.5, {3}'); """, quoted_node_names=True)
        self.pt5 = TreePattern(""" ( '@.name == "A"', '@.dist==0.5{3}'); """, quoted_node_names=True)
        self.pt6 = TreePattern(""" ( '@.name == "C"', '@.dist==0.5{3}'); """, quoted_node_names=True)
        self.pt7 = TreePattern(""" ( '@.name == "A"', '@.dist==0.5{4}'); """, quoted_node_names=True)
        self.pt8 = TreePattern(""" ( '@.name == "A"', '@.dist<0, +'); """, quoted_node_names=True)
        self.pt9 = TreePattern(""" ( '@.name == "A"', '@.dist>0, +'); """, quoted_node_names=True)

        self.pt1_match = True
        self.pt2_match = False
        self.pt3_match = False
        self.pt4_match = True
        self.pt5_match = True
        self.pt6_match = True
        self.pt7_match = False
        self.pt8_match = False
        self.pt9_match = True

    def test_simple_match(self):
        result = len(list(self.pt1.find_match(self.tree))) > 0
        self.assertEqual(result, self.pt1_match)

    def test_simple_match_false_property(self):
        result = len(list(self.pt2.find_match(self.tree))) > 0
        self.assertEqual(result, self.pt2_match)

    def test_relax_with_property(self):
        result = len(list(self.pt3.find_match(self.tree))) > 0
        self.assertEqual(result, self.pt3_match)

    def test_relax_strict_number_with_property_1(self):
        result = len(list(self.pt4.find_match(self.tree))) > 0
        self.assertEqual(result, self.pt4_match)

    def test_relax_strict_number_with_property_2(self):
        result = len(list(self.pt5.find_match(self.tree))) > 0
        self.assertEqual(result, self.pt5_match)

    def test_relax_strict_number_with_property_3(self):
        result = len(list(self.pt6.find_match(self.tree))) > 0
        self.assertEqual(result, self.pt6_match)

    def test_relax_strict_number_with_property_no_match(self):
        result = len(list(self.pt7.find_match(self.tree))) > 0
        self.assertEqual(result, self.pt7_match)

    def test_relax_with_property_no_match(self):
        result = len(list(self.pt8.find_match(self.tree))) > 0
        self.assertEqual(result, self.pt8_match)

    def test_relax_with_property_all_match(self):
        result = len(list(self.pt9.find_match(self.tree))) > 0
        self.assertEqual(result, self.pt9_match)


if __name__ == '__main__':
    unittest.main()
