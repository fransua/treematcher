from itertools import permutations
from ete3 import NCBITaxa
import six
import re
import ast

from ete3 import PhyloTree, Tree, NCBITaxa

class TreePattern(Tree):

    _syntax_tuples = [
        ("@", "__target"),
        ("__target.contains_species", "[n.species for n in temp_leaf_cache[__target]]"),
        ("__target.leaves",  "[n.name for n in temp_leaf_cache[__target]]"),
        ("__target.size", "len(temp_leaf_cache[__target])"),

    ]

    def __init__(self, *args, **kargs):
        kargs["format"] = 1

        self.cache_flag = 0

        is_exact = False

        if len(args) > 0:
            pattern_string = args[0].strip()

            # search exactly as Newick tree, no additional pattern
            if pattern_string.find('@') == -1:
                is_exact = True
            else:
                leaf_searches = ['@.leaves', "@.contains_species", "@.contains", "@.size"]

                for search_type in leaf_searches:
                    if pattern_string.find(search_type) != -1:
                        self.cache_flag = 1

        Tree.__init__(self, *args, **kargs)
        self.temp_leaf_cache={}

        for n in self.traverse():
            if n.name != "NoName":
                self._parse_constraint(n, is_exact)
            else:
                n.constraint = None


    def constrain_match(self, __target, local_vars=None):

        if not local_vars:
            local_vars = {}
        local_vars.update({"__target": __target, "self": __target,
                           "temp_leaf_cache": self.temp_leaf_cache,
                           "smart_lineage": self.smart_lineage,
                           "ncbi": NCBITaxa()})

        try:
            st = eval(self.constraint, local_vars) if self.constraint else True  # eval string as python code

            st = bool(st)  # note that bool of any string returns true
        except ValueError:
                raise ValueError("The following constraint expression did not return boolean result: %s BUT %s" %
                                 (self.constraint, st))
        except (AttributeError, IndexError) as err:
                print('Warning: Constraint evaluation failed at ' + str(__target) + ' with error "' + str(err) + '"')
                st = False

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
        """ create cache  """
        self.temp_leaf_cache = tree.get_cached_content()

        # pass by reference, so dictionary is available on all nodes
        for node in self.traverse():
            node.temp_leaf_cache = self.temp_leaf_cache

        # notes on how to add a custom cache, will do for distance2node
        # self.add_feature("leaves", [node for node in self.get_leaf_names()])
        # temp_cache=tree.get_cached_content(store_attr="leaves", leaves_only=False)
        # this is equivalent to t.get_cached_content(store_attr="name"), just using it here as an example

    def find_match(self, tree, local_vars, maxhits=1):

        if self.cache_flag == 1:
            self.preprocess(tree)

        num_hits = 0
        for node in tree.traverse("preorder"):

            if self.is_match(node, local_vars=local_vars):
                num_hits += 1
                yield node

            if maxhits is not None and num_hits >= maxhits:
                break




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
                    print("Error in syntax dictionary iteration at keyword: " + str(keyword) + "and value: " + python_code)

            if ".lineage" in node.constraint:
                node.constraint = self.smart_lineage(node.constraint)



        return

    def smart_lineage(self, constraint):
        """
        Returns list of lineage tax ids if taxid is searched, otherwise returns names in lineage.
        For example, if a string is given before the "in @.linage" in a query, will get names instead of tax ids.
        Note that names for genus rank and higher ranks must be capitalized.
        """

        #if constraint contains .chilren[0].search, grab that and replace

        parsedPattern = ast.parse(constraint, mode='eval')

        #look in pattern at abstract syntax tree, the left sibling to @.lineage

        lineage_node = [n for n in ast.walk(parsedPattern)
                        if hasattr(n, 'comparators') and type(n.comparators[0]) == ast.Attribute
                        and n.comparators[0].attr == "lineage"]

        index = 0
        for lineage_search in lineage_node:
            if hasattr(lineage_node[index].left,'s'):
                # use re to find  and retrieve what is between __target and .lineage
                found_target = (re.search(r'__target[^ ]*\.lineage', constraint).span())
                extracted_target = constraint[found_target[0] : found_target[1]]

                syntax = "(ncbi.get_taxid_translator(" + \
                         str(extracted_target) + ")).values()"
                if index == 0:
                    constraint = constraint.replace(str(extracted_target), syntax, 1)
                else:
                    # constraint = re.sub(r'^((.*?__target\.lineage.*?){{{index}}})__target\.lineage'.
                    #                    format(index=index),
                    #                    r'\1' + syntax, constraint)

                    # constraint = re.sub(r'^((.*?' + extracted_target + r'.*?){{{index}}})' + extracted_target +
                    #                 .format(index=index),
                    #             r'\1' + syntax, constraint)

                    constraint = re.sub(r'^((.*?' + extracted_target + r'.*?){' + str(index) + r'})' + extracted_target,
                             r'\1' + syntax, constraint)

            # else lineage_node[0].left.n exists and the prefix is a number, don't change anything

            index += 1



        return constraint


def test():
    pattern0 = """ ('   @.dist == 1 and "Gallus_gallus_1" in @.leaves ');"""
    pattern1 = """( '  @.dist >= 0.5 ' , ' @.dist<2  ')
       '    "Pan_troglodytes_1" in @.leaves and "Homo_sapiens_1" in @.children[0] or
       "Pan_troglodytes_1" in @.leaves and "Homo_sapiens_1" in @.children[1]';"""

    pattern0 = TreePattern(pattern0, format=8, quoted_node_names=True)
    pattern1 = TreePattern(pattern1, format=8, quoted_node_names=True)

    tree = PhyloTree(
        "((((Anolis_carolinensis_1:1, Gallus_gallus_1:1), (Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_1:1)primates)primates), ((Danio_rerio_1:1, (Xenopus_laevis_1:1, Anolis_carolinensis_1:1)), Saccharomyces_cerevisiae_2:1)), Saccharomyces_cerevisiae_1:1);",
        format=1)

    print(list(pattern0.find_match(tree, None, maxhits=None)))
    print(len(list(pattern0.find_match(tree, None, maxhits=None))))
    print(len(list(pattern0.find_match(tree, None, maxhits=1))))
    print(len(list(pattern0.find_match(tree, None))))


    print(list(pattern1.find_match(tree, None)))
    print(len(list(pattern1.find_match(tree, None))))


if __name__ == "__main__":
    test()
