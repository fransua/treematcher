from itertools import permutations
from ete3 import NCBITaxa
import six
import re
import ast

from ete3 import PhyloTree, Tree, NCBITaxa


class TreePattern(Tree):
    """
        class for presenting a tree of constrained patterns. A way to represent patterns to search for within trees.

        :param args: Pattern to be searched
        :param kargs: Arguments passed to Tree Class

    """
    _syntax_tuples = [
        ("@", "__target"),
        ("__target.leaves",  "get_cached_attr('name', __target)"),
        ("__target.size", "len(get_cached_attr('name', __target))"),
        #("__target.contains_species", "get_cached_attr('species', __target)"),
        #("__target.size", "len(temp_leaf_cache[__target])"),
        # ("__target.contains_species", "[n.species for n in temp_leaf_cache[__target]]"),
        #("__target.leaves", "[n.name for n in temp_leaf_cache[__target]]"),

    ]

    def __init__(self, *args, **kargs):

        kargs["format"] = 1

        self.cache_flag = 0

        is_exact = False

        self.all_node_cache = None

        try:
            self._extra_functions = kargs.pop('functions')

        except KeyError:
            self._extra_functions = None

        if len(args) > 0:
            pattern_string = args[0].strip()

            if pattern_string.find('@') == -1:
                is_exact = True
            else:
                searches = ['.leaves', ".contains_species", ".contains", ".size"]

                for search_type in searches:
                    if pattern_string.find(search_type) != -1:
                        self.cache_flag = 1

        #Tree.__init__(self, *args, **kargs)
        super(TreePattern, self).__init__(*args, **kargs)
        self.temp_leaf_cache = {}
        # all of this only happens in the root node
        if len(args)>0:

            for n in self.traverse():
                n._extra_functions = self._extra_functions

                if n.name != "NoName":
                    self._parse_constraint(n, is_exact)
                else:
                    n.constraint = None


    def constrain_match(self, __target, local_vars=None):
        """
        Evaluate constraint on a single node. Checks whether a single node matches a single constraint.
        :param __target: represents a node you are looking to target. Replaces @ in a query.
        :param local_vars: Dictionary of treematcher class variables and functions for constraint evaluation
        :return: returns a boolean value of True if a match is found, otherwise False.
        """
        if not local_vars:
            local_vars = {}
        local_vars.update({"__target": __target, "self": __target,
                           "temp_leaf_cache": self.temp_leaf_cache,
                           "smart_lineage": self.smart_lineage,
                           "ncbi": NCBITaxa(),
                           "pattern": self,
                           "get_cached_attr": self.get_cached_attr
                           })
        if self._extra_functions:
            local_vars.update(self._extra_functions)

        try:
            st = eval(self.constraint, local_vars) if self.constraint else True

            st = bool(st)  # note that bool of any string returns true
        except ValueError:
            raise ValueError("The following constraint expression did not return boolean result: %s BUT %s" %
                             (self.constraint, st))
        except (AttributeError, IndexError) as err:
            print('Warning: Constraint evaluation failed at ' + str(__target) + ' with error "' + str(err) + '"')
            st = False

        return st


    def is_match(self, node, local_vars=None):

        """
        Check all constraints on a tree. Permutes the tree and checks that all constrinats return True. CHecks whther a single node and its children match a given pattern
        :param node: A tree (root node) to be searched for a given pattern
        :param local_vars:  Dictionary of treematcher class variables and functions for constraint evaluation
        :return: True is a match has been found, otherwise False
        """
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
        """
        :param tree: Patern tree. A cache is made available on all nodes of the tree when searches require
                    repetitive traversal.

        """
        self.temp_leaf_cache = tree.get_cached_content()

        self.all_node_cache = tree.get_cached_content(leaves_only=False)

        # pass by reference, so dictionary is available on all nodes
        for node in self.traverse():
            node.temp_leaf_cache = self.temp_leaf_cache
            node.all_node_cache = self.all_node_cache

    def find_match(self, tree, local_vars, maxhits=1):

        """

        :param tree: tree to be searched for a matching pattern.
        :param local_vars:  Dictionary of treematcher class variables and functions for constraint evaluation
        :param maxhits: Number of matches to be searched for.
        :param None maxhits: Pattern search loop will continue until all matches are found.
        """
        if self.cache_flag == 1:
            self.preprocess(tree)

        num_hits = 0
        num = 0
        for node in tree.traverse("preorder"):
            num += 1

            if self.is_match(node, local_vars=local_vars):
                print num, '- YES'
                print node.write(format=9)
                # print node
                num_hits += 1
                yield node
            else:
                print num, '- NOP'
                print node.write(format=9)
                # print node
            print '-'*60
            if maxhits is not None and num_hits >= maxhits:
                break


    def _parse_constraint(self, node, is_exact=False):
        """
        Function for replacing keywords with python code in a query pattern.

        :param node: The pattern be searched for in the query.
        :param is_exact: Designates if pattern is in standard Newick format or contains treematcher syntax.
        :param True is_exact: Pattern is a tree standard Newick format, search for exact match of tree.
        :param False is_exact:  Pattern contains treematcher syntax and requires keyword to python conversion.
        """
        node.constraint = node.name

        #builtins = set(PhyloTree.__dict__.keys() + Tree.__dict__.keys())
        

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

            #for attrib in re.findall("__target.([A-Za-z_0-9]+)\(", node.constraint):
            #    if attrib not in builtins:
            #        node.constraint = re.sub("__target.([A-Za-z_0-9]+)\(([^)]+)",
            #                                 "\\1(__target, pattern, (\\2)", node.constraint)

            #print "self._extra_functions are:", self._extra_functions
            #Search for custom functions
            for custom_function in self._extra_functions.keys():
                #node.constraint = re.findall(custom_function + "([A-Za-z_0-9]+)\(", node.constraint)
                node.constraint = re.sub(custom_function + "\(__target",
                                         custom_function + "(__target, pattern", node.constraint)

            #custom function of stuff tpat=whatever

            if ".lineage" in node.constraint:
                node.constraint = self.smart_lineage(node.constraint)

            # if "contains_species" in node.constraint:
            #     node.constraint = self.contains_species(node.constraint, node)

    def get_cached_attr(self, attr_name, node):
        """
        :param attr_name: any attribute cached in tree nodes (e.g., species, name, dist, support, etc.)
        :param node: The pattern tre containing the cache
        :return: cached values for the requested attribute (e.g., Homo sapiens, Human, 1.0, etc.)
        """
        values = list(getattr(n, attr_name, None) for n in self.all_node_cache[node])
        #values.discard(None)
        return values

    def smart_lineage(self, constraint):

        """
        For example, if a string is given before the "in @.linage" in a query, will get names instead of tax ids.
        Note that names for genus rank and higher ranks must be capitalized. Function should work for
         constraint that contains something besides the given target node  (e.g., @.chilren[0].lineage)
        :param constraint: The entire pattern being searched for which includes @.lineage.
        :return:  Returns list of lineage tax ids if taxid is searched, otherwise returns names in lineage.

        """
        parsedPattern = ast.parse(constraint, mode='eval')

        lineage_node = [n for n in ast.walk(parsedPattern)
                        if hasattr(n, 'comparators') and type(n.comparators[0]) == ast.Attribute
                        and n.comparators[0].attr == "lineage"]

        index = 0
        for lineage_search in lineage_node:
            if hasattr(lineage_node[index].left,'s'):
                # retrieve what is between __target and .lineage
                found_target = (re.search(r'__target[^ ]*\.lineage', constraint).span())
                extracted_target = constraint[found_target[0]: found_target[1]]

                syntax = "(ncbi.get_taxid_translator(" + \
                         str(extracted_target) + ")).values()"
                if index == 0:
                    constraint = constraint.replace(str(extracted_target), syntax, 1)
                else:

                    constraint = re.sub(r'^((.*?' + extracted_target + r'.*?){' + str(index) + r'})' + extracted_target,
                             r'\1' + syntax, constraint)

            index += 1



        return constraint


def contains_species(__target, tpat, constraint):
    tpat.preprocess(__target)
    if isinstance(constraint, six.string_types):
        species_list = [constraint]
    else:
        species_list=constraint
    cached_species = tpat.get_cached_attr('species', __target)
    result = all(sp in cached_species for sp in species_list)

    return result


def contains(__target, tpat, constraint):
    tpat.preprocess(__target)
    if isinstance(constraint, six.string_types):
        name_list = [constraint]
    else:
        name_list = constraint
    cached_name = tpat.get_cached_attr('name', __target)
    result = all(name in cached_name for name in name_list)
    return result


def number_of_species(__target, tpat, constraint):
    tpat.preprocess(__target)
    cached_species = tpat.get_cached_attr('species', __target)
    result = len(cached_species)
    print("result is", result)

    return result


def num_of_leaves(__target):
    tpat.preprocess(__target)
    cached_name = tpat.get_cached_attr('name', __target)
    result = len(cached_name)
    print("result is", result)

    return result

def is_duplication(__target):
    # checks node.evol_type inferred with PhyloTree.get_descendant_evol_events()
    pass

def is_speciation(__target):
    # checks
    pass


def has_match():
    # call find match and return true if the target tree has at least one match
    pass



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

def test2():
    t = PhyloTree(
        "((((Human_1, Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1);")
    t.set_species_naming_function(lambda node: node.name.split("_")[0])
    #ntrees, ndups, sptrees = t.get_speciation_trees()
    #print "Found %d species trees and %d duplication nodes" % (ntrees, ndups)
    #for spt in sptrees:
    #    print spt
    #print t
    #pattern = """( ' "Human" in @.contains_species and "Chimp" in @.contains_species   '); """
    #pattern = TreePattern(pattern, format=8, quoted_node_names=True)
    #print(len(list(pattern.find_match(t, None, maxhits=None))))

    pattern1 = """( 'contains(@, ("Chimp_2", "Chimp_3"))',  'number_of_species(@, 1) '); """
    tp1 = TreePattern(pattern1, format=8, quoted_node_names=True,
                      functions={'contains': contains,
                                 "contains_species": contains_species,
                                 "number_of_species": number_of_species})
    print(len(list(tp1.find_match(t, None, maxhits=None))))


if __name__ == "__main__":
    #test()
    test2()
