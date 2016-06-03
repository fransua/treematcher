from itertools import permutations
from string import strip
import re
import sys

from ete3 import PhyloTree, Tree

class TreePattern(Tree):

    _syntax_tuples = [
        ("@", "__target"),
        ("__target.leaves", "[node.name for node in __target.get_leaves()]"),
        ("__target.children", "[n.name for n in __target.children]"),
        ("is greater than", ">"),
        ("greater than", ">"),
        (" or equal to", "="),
        ("is less than", "<"),
        ("less than", "<"),
        ("Children", "[n.name for n in __target.children]"),
        ("Distance", "__target.dist"),
        ("Leaves", "[node.name for node in __target.get_leaves()]"),
        ("Species", "__target.species"),
        (" is ", " == "),
        ]

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

        node.constraint = node.name
        # use regular expressions turn multiple spaces to single space
        node.constraint = re.sub("\s+", " ", node.constraint)

        for keyword, python_code in self._syntax_tuples:
            try:
                node.constraint = node.constraint.replace(keyword, python_code)
            except (KeyError, ValueError):
                print "Error in syntax dictionary iteration at keyword: " + str(keyword) + "and value: " + python_code

        return

def length(txt):
    return len(txt)

def test_basic():
    """
        tests basic ete parameters are working
    """
    pattern = """
            (
            'len(@.children) > 2 and @.name in ("hello","bye") '
            )
            '(len(@.name) < 3) and @.dist >= 0.5';
            """

    pattern = TreePattern(pattern, format=8, quoted_node_names=True)

    print pattern

    tree = Tree("(hello,(1,2,3)kk)pasa:1;", format=1)
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree?:", pattern.find_match(tree, None)

    tree = Tree("((kk,(1,2,3)bye)y:1, NODE);", format=1)
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree?:", pattern.find_match(tree, None)

    tree = Tree("(((1,2,3)bye)y:1, NoName);", format=1)
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree?:", pattern.find_match(tree, None)


def test_syntax():
    """
        tests syntax-to-python conversion is working
    """
    # " Saccharomyces_cerevisiae_1" in Children
    # "Homo_sapiens_1" in @.children

    pattern1 = """
        ( '  @.dist >= 0.5 ' , ' @.dist<2  ')
        '    "Pan_troglodytes_1" in @.leaves and "Homo_sapiens_1" in @.children '
        ;
        """

    pattern2 = """
        ( '  Distance greater than or equal to 0.5 ' , ' Distance less than 2 ' )
        '     "Pan_troglodytes_1" in Leaves and "Homo_sapiens_1" in Children'
        ;
        """

    pattern1 = TreePattern(pattern1, format=8, quoted_node_names=True)
    pattern2 = TreePattern(pattern2, format=8, quoted_node_names=True)

    print pattern1
    print pattern2


    tree = PhyloTree("((((Anolis_carolinensis_1:1, Gallus_gallus_1:1), (Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_1:1))), ((Danio_rerio_1:1, (Xenopus_laevis_1:1, Anolis_carolinensis_1:1)), Saccharomyces_cerevisiae_2:1)), Saccharomyces_cerevisiae_1:1);", format=1)
    tree.set_species_naming_function(lambda n: n.name.split("_")[1] if "_" in n.name else '')
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree?:", pattern1.find_match(tree, None)
    print "Pattern without symbols matches tree?:", pattern2.find_match(tree, None)


    tree = PhyloTree("((((Anolis_carolinensis_1:1, Gallus_gallus_1:1), (Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_2:1))), ((Danio_rerio_1:1, (Xenopus_laevis_1:1, Anolis_carolinensis_1:1)), Saccharomyces_cerevisiae_2:1)), Saccharomyces_cerevisiae_1:1);", format=1)
    tree.set_species_naming_function(lambda n: n.name.split("_")[1] if "_" in n.name else '')
    #print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree missing leaf?:", pattern1.find_match(tree, None)
    print "Pattern without symbols matches tree missing leaf?:", pattern2.find_match(tree, None)

def test_species():
    """
        tests if @.species is working
    """

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
    print "Pattern without symbols matches tree?:", pattern2.find_match(tree, None)

    tree = PhyloTree("((((Anolis_carolinensis_1:1, Gallus_gallus_1:1), (Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_1:1))), ((Danio_rerio_1:1, (Xenopus_laevis_1:1, Anolis_carolinensis_1:1)), Saccharomyces_cerevisiae_2:1)), Saccharomyces_cerevisiae_1:1);", format=1)
    tree.set_species_naming_function(lambda node: (node.name.split("_")[0] + " " + node.name.split("_")[1]) if "_" in node.name else '')
    #print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches wrong species tree?:", pattern1.find_match(tree, None)
    print "Pattern without symbols matches wrong species tree?:", pattern2.find_match(tree, None)

    tree = PhyloTree("((((Anolis_carolinensis_1:1, Gallus_gallus_1:1), (Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_2:1))), ((Danio_rerio_1:1, (Xenopus_laevis_1:1, Anolis_carolinensis_1:1)), Saccharomyces_cerevisiae_2:1)), Saccharomyces_cerevisiae_1:1);", format=1)
    tree.set_species_naming_function(lambda n: n.name.split("_")[1] if "_" in n.name else '')
    #print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree missing leaf?:", pattern1.find_match(tree, None)
    print "Pattern without symbols matches tree missing leaf?:", pattern2.find_match(tree, None)

def test_custom_functions():
    """
        tests all custom functions are working
    """
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



if __name__ == "__main__":

    test_basic()
    test_syntax()
    test_species()
    #test_custom_functions()
    #test_evol_events()



################################
########## NOTES ###############

    ##### Potential Issues ######

    # 1) @.species is "sapiens" or "pygmaeus"
    #    is the same as
    #    @.species == "sapiens"or node.name=="pygmaeus"
    #    which may not be what people expect
    # 2) @.species will fail if not all nodes have species
    #
    # 3) if someone were to name their leaves with keywords, they would be modified
    #   need to replace inner strings with temporary variables during the replace


    ##### To Do ######
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
    # may want to store syntax metadata in an external file for easy access and minimal loading required
    #
    # when parsing constraints
    # may want to pull out everything in quotes first or extract some python keywords (in, is, or, and) first?
    #
    # add regular expressions to each keyword as needed
