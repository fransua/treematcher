from itertools import permutations
from string import strip
import re
import sys

from ete3 import PhyloTree, Tree

class TreePattern(Tree):
    def __init__(self, *args, **kargs):
        kargs["format"] = 1
        Tree.__init__(self, *args, **kargs)
        for n in self.traverse():
            if n.name != "NoName":
                # n.constraint = n.name.replace("@", "__target")
                self._parse_constraint(n)
            else:
                n.constraint = None

    def constrain_match(self, __target, local_vars=None):
        # if not self.constraint:
        #   return True
        # redundant true in line 26

        if not local_vars:
            local_vars = {}
        local_vars.update({"__target": __target, "self": __target})
        try:
            st = eval(self.constraint, local_vars) if self.constraint else True  # eval string as python code
            #print __target

            st = bool(st)  # note that bool of any string returns true
        except ValueError: 
                raise ValueError("The following constraint expression did not return boolean result: %s BUT %s" %
                                 (self.constraint, st))

        return st
    
    def is_match(self, node, local_vars=None):
        # Check expected features
        status = self.constrain_match(node, local_vars)
        if status and self.children:
            #print "has children"
            if len(node.children) >= len(self.children):
                # Check all possible comparison between pattern children and
                # and tree node children.
                for candidate in permutations(node.children):
                    sub_status = True
                    for i in range(len(self.children)):
                        st = self.children[i].is_match(candidate[i], local_vars)
                        sub_status &= st
                    status = sub_status
                    if status:
                        break
            else:
                status = False
        return status
    
    def __str__(self):
        return self.get_ascii(show_internal=True, attributes=["constraint"])

    def find_match(self, tree, local_vars):
        for node in tree.traverse("preorder"):
            if self.is_match(node, local_vars=local_vars):
                return True, node
        return False, None

    def _parse_constraint(self,node):
        node.constraint = node.name.replace("@", "__target")\
            .replace("Leaves", "[node.name for node in __target.get_leaves()]")\
            .replace("__target.leaves", "[node.name for node in __target.get_leaves()]")\
            .replace("Distance", "__target.dist")\
            .replace("Species", "__target.species")\
            .replace("less than", "<") \
            .replace(" or equal to", "=") \
            .replace("greater than", ">") \
            .replace(" is ", " == ")

        return


def length(txt):
    return len(txt)


def test_basic():
    custom_functions = {"length":length}


    pattern = """
        (
        'len(@.children) > 2 and @.name in ("hello","bye") '
        )
        '(length(@.name) < 3) and @.dist >= 0.5';
        """

    pattern = TreePattern(pattern, format=8, quoted_node_names=True)

    print pattern

    tree = Tree("(hello,(1,2,3)kk)pasa:1;", format=1)
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree?:", pattern.find_match(tree, custom_functions)

    tree = Tree("((kk,(1,2,3)bye)y:1, NODE);", format=1)
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree?:", pattern.find_match(tree, custom_functions)

    tree = Tree("(((1,2,3)bye)y:1, NODE);", format=1)
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree?:", pattern.find_match(tree, custom_functions)

    tree = Tree("(((1,2,3)bye,kk)y:1, NODE);", format=1)
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree?:", pattern.find_match(tree, custom_functions)


def test_syntax():

    pattern1 = """
        (
        (
        ( '  @.dist >= 0.5 ' , ' @.species in ("sapiens","pygmaeus")  ')
        )
        ' "Pan_troglodytes_1" in @.leaves and @.dist<2 '
        )
        ;
        """

    pattern2 = """
        (
        (
        ( '  Distance greater than or equal to 0.5 ' , ' Species is "sapiens" or Species is "pygmaeus" ')
        )
        ' "Pan_troglodytes_1" in Leaves and Distance less than 2'
        )
        ;
        """

    pattern1 = TreePattern(pattern1, format=8, quoted_node_names=True)
    pattern2 = TreePattern(pattern2, format=8, quoted_node_names=True)

    print pattern1
    print pattern2


    tree = PhyloTree("((((Anolis_carolinensis_1:1, Gallus_gallus_1:1), (Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_1:1))), ((Danio_rerio_1:1, (Xenopus_laevis_1:1, Anolis_carolinensis_1:1)), Saccharomyces_cerevisiae_2:1)), Saccharomyces_cerevisiae_1:1);", format=1)
    tree.set_species_naming_function(lambda n: n.name.split("_")[1] if "_" in n.name else '')
    print tree.get_ascii(attributes=["species", "dist"])
    print "Pattern matches tree?:", pattern1.find_match(tree, None)
    print "Pattern without @ matches tree?:", pattern2.find_match(tree, None)

    tree = PhyloTree("((((Anolis_carolinensis_1:1, Gallus_gallus_1:1), (Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_1:1))), ((Danio_rerio_1:1, (Xenopus_laevis_1:1, Anolis_carolinensis_1:1)), Saccharomyces_cerevisiae_2:1)), Saccharomyces_cerevisiae_1:1);", format=1)
    tree.set_species_naming_function(lambda node: (node.name.split("_")[0] + " " + node.name.split("_")[1]) if "_" in node.name else '')
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches wrong species tree?:", pattern1.find_match(tree, None)
    print "Pattern without @ matches wrong species tree?:", pattern2.find_match(tree, None)

    tree = PhyloTree("((((Anolis_carolinensis_1:1, Gallus_gallus_1:1), (Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_2:1))), ((Danio_rerio_1:1, (Xenopus_laevis_1:1, Anolis_carolinensis_1:1)), Saccharomyces_cerevisiae_2:1)), Saccharomyces_cerevisiae_1:1);", format=1)
    tree.set_species_naming_function(lambda n: n.name.split("_")[1] if "_" in n.name else '')
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree missing leaf?:", pattern1.find_match(tree, None)
    print "Pattern without @ matches tree missing leaf?:", pattern2.find_match(tree, None)


if __name__ == "__main__":
    test_basic()
    test_syntax()


################################
########## NOTES ###############

    ##### Potential Issues ######

    #  Genus is "Homo" or "Pan"
    # is the same as
    # Genus is "Homo" or node.name=="Pan"

    # @.species will fail if not all nodes have species

    ##### To Do ######
    #
    # parse strings using regular expressions (perhaps use an ordered dictionary on keywords?)
    #
    # use keyword Subtree to apply to all nodes children
    #
    # .replace("@.lineage","__target.get_taxonomy")
    # .replace("Lineage", "[(ncbi.get_taxid_translator(lineage))[taxid] for taxid in ncbi.get_lineage(__target.taxid)]")\
    # .replace("Genus","__target.get_taxonomy("Genus")")
    # .replace("Duplication", '__target.evol_event=="D"')\
    #
    # Examples:
    #   "Hominidae" in @.lineage
    #   "Hominidae" in Lineage
    #   Lineage contains("Hominidae") # would require a contains keyword to be defined with re
    #
    #