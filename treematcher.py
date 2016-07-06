from itertools import permutations
import re

from ete3 import PhyloTree, Tree, NCBITaxa

class TreePattern(Tree):

    _syntax_tuples = [
        ("@", "__target"),
        ("__target.size", "len(self.get_tree_root().temp_cache_name[__target])"),
        ("__target.leaves", "[node.name for node in __target.get_leaves()]"),
        ("__target.children", "[n.name for n in __target.children]"),
        ("__target.contains", "self.get_tree_root().temp_cache_name[__target]"),
        #("__target.leaves", "__target.temp_cache_leaves[__target]"),
        #("__target.children", "__target.temp_cache_children[__target]"),
    ]

    def __init__(self, *args, **kargs):
        kargs["format"] = 1
        #initialize cache flags
        self.cache_flag_dict = {
            "@.contains": False,
            "@.leaves": False,
            "@.children": False,
            #"@.contains_species": False,
        }

        is_exact = False

        if len(args) > 0:
            pattern_string = args[0].strip()

            # search exactly as Newick tree, no additional pattern
            if pattern_string.find('@') == -1:
                is_exact = True
            else:
             searches = ['@.contains', '@.leaves', '@.children']
             for search_type in searches:
                if pattern_string.find(search_type) != -1:
                    self.cache_flag_dict[search_type] = True

        Tree.__init__(self, *args, **kargs)


        for n in self.traverse():

            if n.name != "NoName":
                self._parse_constraint(n, is_exact)
            else:
                n.constraint = None


    def constrain_match(self, __target, local_vars=None):

        if not local_vars:
            local_vars = {}
        local_vars.update({"__target": __target, "self": __target})

        try:
            st = eval(self.constraint, local_vars) if self.constraint else True  # eval string as python code

            st = bool(st)  # note that bool of any string returns true
        except ValueError:
                raise ValueError("The following constraint expression did not return boolean result: %s BUT %s" %
                                 (self.constraint, st))

        return st


    def is_match(self, node, local_vars=None):
        # Check expected features

        status = self.constrain_match(node, local_vars)

        if status and self.children:
            if len(node.children) >= len(self.children):
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


    def preprocess(self, tree):
        root = tree.get_tree_root()
        if self.cache_flag_dict['@.contains']:
            root.add_feature("temp_cache_name", tree.get_cached_content(store_attr="name"))

        #if self.cache_flag_dict['@.leaves']:
        #    for n in tree:
        #        n.add_feature("temp_cache_leaves", [node.name for node in n.get_leaves()])

        #if self.cache_flag_dict['@.children']:
        #    for n in tree:
        #        n.add_feature("temp_cache_children", [node.name for node in n.get_children()])



    def find_match(self, tree, local_vars):
        self.preprocess(tree)

        for node in tree.traverse("preorder"):
            if self.is_match(node, local_vars=local_vars):
                return True, node
        return False, None

    def _parse_constraint(self, node, is_exact=False):
        node.constraint = node.name

        if is_exact and node.name != '':
            node.constraint = "__target.name==" + "'" + str(node.constraint) + "'"
        else:
            node.constraint = node.name
            # turn multiple spaces to single space

            node.constraint = re.sub("\s+", " ", node.constraint)

            for keyword, python_code in self._syntax_tuples:
                try:
                    node.constraint = node.constraint.replace(keyword, python_code)
                except (KeyError, ValueError):
                    print "Error in syntax dictionary iteration at keyword: " + str(keyword) + "and value: " + python_code

        return


def test():
    pattern0 = """ ('   @.size == 1 and "Gallus_gallus_1" in @.contains ');"""
    pattern1 = """( '  @.dist >= 0.5 ' , ' @.dist<2  ')
        '    "Pan_troglodytes_1" in @.leaves and "Homo_sapiens_1" in @.children ';"""

    pattern0 = TreePattern(pattern0, format=8, quoted_node_names=True)
    pattern1 = TreePattern(pattern1, format=8, quoted_node_names=True)

    tree = Tree(
        "((((Anolis_carolinensis_1:1, Gallus_gallus_1:1), (Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_1:1)primates)primates), ((Danio_rerio_1:1, (Xenopus_laevis_1:1, Anolis_carolinensis_1:1)), Saccharomyces_cerevisiae_2:1)), Saccharomyces_cerevisiae_1:1);",
        format=1)
    #print tree
    print pattern0.find_match(tree, None)
    print pattern1.find_match(tree, None)[1].name

if __name__ == "__main__":
    test()
