from itertools import permutations
from ete3 import NCBITaxa
import six
import re
import ast

from ete3 import PhyloTree, Tree, NCBITaxa


class TreePattern(Tree):
    """
        Class for presenting a tree of constrained patterns. A way to represent patterns to search for within trees.
        :param args: Pattern to be searched.
        :param kargs: Custom functions and arguments passed to Tree Class.
    """
    _syntax_tuples = [
        ("@", "__target"),
        ("__target.size", "len(temp_leaf_cache[__target])"),
        ("__target.contains_species", "[n.species for n in temp_leaf_cache[__target]]"),
        ("__target.leaves", "[n.name for n in temp_leaf_cache[__target]]"),
        #("__target.leaves", "get_cached_attr('name', __target)"), # no longer returns leaves, need similar function
        #("__target.size", "len(get_cached_attr('name', __target))"),
        #("__target.contains_species", "get_cached_attr('species', __target)"),
    ]

    def __init__(self, *args, **kargs):
        kargs["format"] = 1

        is_exact = False

        self.cache_flag = 0
        self.all_node_cache = None

        self.evol_events_flag = 0
        self.evol_events = None
        self.shortcut_functions = []



        try:
            self._extra_functions = kargs.pop('functions')

        except KeyError:
            self._extra_functions = None


        if len(args) > 0:
            pattern_string = args[0].strip()

            if pattern_string.find('@') == -1:
                is_exact = True
            else:
                cached_searches = ['.leaves', ".contains_species", ".contains", ".size"]

                shortcut_functions=["contains_species", "contains_duplication", "contains_speciation", "contains"]



                for search_type in cached_searches:
                    if pattern_string.find(search_type) != -1:
                        self.cache_flag = 1
                for function in shortcut_functions:
                    if pattern_string.find(function + '(') != -1:
                        self.shortcut_functions += [function]


                if pattern_string.find("contains_duplication")!= -1 or pattern_string.find("contains_speciation") != -1:
                    self.evol_events_flag = 1


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


    def constrain_match(self, target_node, local_vars=None):
        """
        Evaluate constraint on a single node.
        :param target_node: Node evaluated
        :param local_vars: Dictionary of treematcher class variables and functions for constraint evaluation
        :return: Returns True if a match is found, otherwise False.
        """
        if not local_vars:
            local_vars = {}
        local_vars.update({"__target": target_node, "self": target_node,
                           "temp_leaf_cache": self.temp_leaf_cache,
                           "smart_lineage": self.smart_lineage,
                           "ncbi": NCBITaxa(),
                           "pattern": self,
                           "get_cached_attr": self.get_cached_attr,
                           "contains_species": contains_species,
                           "contains_duplication": contains_duplication,
                           "contains_speciation": contains_speciation,
                           "contains": contains,
                           "number_of_species": number_of_species,
                           "number_of_leaves": number_of_leaves
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
            print('Warning: Constraint evaluation failed at ' + str(target_node) + ' with error "' + str(err) + '"')
            st = False

        return st


    def is_match(self, node, local_vars=None):
        """
        Check all constraints on the subtree under node.
        :param node: A tree (node) to be searched for a given pattern.
        :param local_vars:  Dictionary of TreePattern class variables and functions for constraint evaluation.
        :return: True if a match has been found, otherwise False.
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
        """
        Overrides method in TreeNode
        """
        return self.get_ascii(show_internal=True, attributes=["constraint"])

    def preprocess(self, tree):
        """
        Creates cache for attributes that require multiple tree traversal and makes the cache available on all nodes.
        :param tree: Pattern to search for.
        """
        self.temp_leaf_cache = tree.get_cached_content()
        #self.evol_events = tree.get_descendant_evol_events()

        self.all_node_cache = tree.get_cached_content(leaves_only=False)

        # pass by reference, so dictionary is available on all nodes
        for node in self.traverse():
            node.temp_leaf_cache = self.temp_leaf_cache
            node.all_node_cache = self.all_node_cache

    def find_match(self, tree, local_vars, maxhits=1):
        """
        A pattern search continues until the number of specified matches are found.
        :param tree: tree to be searched for a matching pattern.
        :param local_vars:  Dictionary of TreePattern class variables and functions for constraint evaluation
        :param maxhits: Number of matches to be searched for.
        :param None maxhits: Pattern search will continue until all matches are found.
        """
        if self.cache_flag == 1 or self.evol_events_flag == 1:
            self.preprocess(tree)

        num_hits = 0
        num = 0
        for node in tree.traverse("preorder"):
            num += 1

            if self.is_match(node, local_vars=local_vars):
                #print(num, '- YES')
                #print(node.write(format=9))
                #print(node)
                num_hits += 1
                yield node
            #else:
                #print(num, '- NOP')
                #print(node.write(format=9))
                #print(node)
            #print('-'*60)
            if maxhits is not None and num_hits >= maxhits:
                break

    def has_match(self, tree, local_vars):
        """
        calls find_match and return true if the target tree has at least one match
        """
        if len(self.find_match(tree, None, maxhits=1))>=1:
            return True
        else:
            return False

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


            if self._extra_functions is not None:
                for custom_function in self._extra_functions.keys():
                    node.constraint = re.sub(custom_function + "\(__target",
                                             custom_function + "(__target, pattern", node.constraint)
            if len(self.shortcut_functions)>0:
                for custom_function in self.shortcut_functions:
                    node.constraint = re.sub(custom_function + "\(__target",
                                             custom_function + "(__target, pattern", node.constraint)

            #custom function of stuff tpat=whatever

            if ".lineage" in node.constraint:
                node.constraint = self.smart_lineage(node.constraint)


    def get_cached_attr(self, attr_name, node):
        """
        Access cached attributes on trees.
        :param attr_name: any attribute cached in tree nodes (e.g., species, name, dist, support, etc.)
        :param node: The pattern tree containing the cache
        :return: cached values for the requested attribute (e.g., Homo sapiens, Human, 1.0, etc.)
        """
        if self.all_node_cache is None:
            self.preprocess(node.get_tree_root())
        values = list(getattr(n, attr_name, None) for n in self.all_node_cache[node])
        #values.discard(None)
        return values


    def smart_lineage(self, constraint):
        """
        Get names instead of tax ids if a string is given before the "in @.linage" in a query.
        Otherwise, returns Taxonomy ids. Function also works for constraint that contains something besides the given target node
        (e.g., @.children[0].lineage).
        :param constraint: Internal use.
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


def contains_species(target_node, pattern, species_list):
    """
    Shortcut function to find the species at a node and any of it's descendants.
    :param target_node: Node to be evaluated, given as @.
    :param pattern: Internal use
    :param species_list: list of species being searched for
    :return: True if all species in list are found, otherwise False.
    """

    if isinstance(species_list, six.string_types):
        species_list = [species_list]

    cached_species = pattern.get_cached_attr('species', target_node)
    result = all(sp in cached_species for sp in species_list)
    return result


def contains(target_node, pattern, name_list):
    """
    Shortcut function to find the name(s) within a node and any of it's descendants.
    :param target_node: Node to be evaluated, given as @.
    :param pattern: Internal use
    :param name_list: List of names to be searched for at a node.
    :return: True if all names in name_list are found, otherwise False.
    """

    if isinstance(name_list, six.string_types):
        name_list = [name_list]

    cached_name = pattern.get_cached_attr('name', target_node)
    result = all(name in cached_name for name in name_list)
    return result


def number_of_species(target_node, pattern, num):
    """
    Shortcut function to find the number of species within a node and any of it's descendants.
    :param target_node: Node to be evaluated, given as @.
    :param pattern: Internal use.
    :param num: number of species searched for at a node.
    :return: True if number of descendant species is equal to num, otherwise False.
    """
    cached_species = pattern.get_cached_attr('species', target_node)
    if num == len(cached_species):
        result = True
    else:
        result = False

    return result


def number_of_leaves(target_node, pattern, num):
    """
    Shortcut function to find the number of leaves within a node and any of it's descendants.
    :param target_node: Node to be evaluated, given as @.
    :param pattern: Internal use.
    :param num: Number of descendant leaves to be searched for.
    :return: True if all number of leaves matches num, otherwise False.
    """
    cached_name = pattern.get_cached_attr('name', target_node)
    if num == len(cached_name):
        result = True
    else:
        result = False

    return result


def contains_duplication(target_node, pattern):
    """
        Shortcut function to find the number of duplication events at or below a node.
        :param target_node: Node to be evaluated, given as @.
        :return: True if node is a duplication, otherwise False.
    """
    events = pattern.get_cached_attr('evoltype', target_node)
    return(events.count('D'))


def contains_speciation(target_node, pattern):
    """
        Shortcut function to find the number of speciation events at or below a node.
        :param target_node: Node to be evaluated, use @
        :return: True if node is a duplication, otherwise False.
    """
    events = pattern.get_cached_attr('evoltype', target_node)
    return(events.count('S'))


def test():
    t = PhyloTree(
        "((((Human_1, Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1);")
    t.set_species_naming_function(lambda node: node.name.split("_")[0])
    t.get_descendant_evol_events()

    pattern0 = """ "Human_2" in @.leaves; """
    pattern1 = """( 'contains(@, ("Chimp_2", "Chimp_3"))'); """

    pattern0 = TreePattern(pattern0, format=8, quoted_node_names=True)
    pattern1 = TreePattern(pattern1, format=8, quoted_node_names=True)


    #should return 5 results
    print(len(list(pattern0.find_match(t, None, maxhits=None))))

    #should return 4 results
    print(len(list(pattern1.find_match(t, None, maxhits=None))))

def test_tutorial():

    tree = Tree("((sample_1,sample_2)ancestor_a,(sample_1,sample_2)ancestor_b)root;", format=8)
    print tree
    border= "\n"+ "#" * 110 + "\n"

    #######################################################
    print(border)
    print("Find a node named sample_1.")
    print(border)
    #######################################################

    pattern1_v1 = """ @.name=="sample_1" ; """
    pattern1_v1 = TreePattern(pattern1_v1)  # create a TreePattern Instance
    solution = pattern1_v1.find_match(tree, None)
    print(pattern1_v1, list(solution))

    pattern1_v2 = """ sample_1 ; """
    pattern1_v2 = TreePattern(pattern1_v2)  # create a TreePattern Instance
    #print("Pattern 1 version 2 is", pattern1_v2)
    solution = pattern1_v2.find_match(tree, None)
    print(pattern1_v2, list(solution))

    pattern1_v3 = """ "sample_1" ; """  # Note that single quotes on the inside will not work here
    pattern1_v3 = TreePattern(pattern1_v3, quoted_node_names=True)  # create a TreePattern Instance
    #print("Pattern 1 version 3 is", pattern1_v3)
    solution = pattern1_v3.find_match(tree, None)
    print(pattern1_v3, list(solution))

    pattern1_v4 = ' "sample_1" ; '
    pattern1_v4 = TreePattern(pattern1_v4, quoted_node_names=True)  # create a TreePattern Instance
    #print("Pattern 1 version 4 is", pattern1_v4)
    solution = pattern1_v4.find_match(tree, None)
    print(pattern1_v4, list(solution))

    # List all pattern 1 matches
    solution = list(pattern1_v2.find_match(tree, None, maxhits=None))
    print("All solutions for pattern 1 are:")
    print(solution)

    # Find the total number of pattern2 matches
    solution = len(list(pattern1_v2.find_match(tree, None, maxhits=None)))
    print("The number of solutions for pattern 1 are:")
    print(solution)

    #######################################################
    print(border)
    print("Find a tree where sample_1 and sample_2 are siblings.")
    print(border)
    #######################################################

    pattern2_v1 = """ @.children[0].name=="sample_1" and @.children[1].name=="sample_2" ; """  # comma is used separate sibling nodes
    pattern2_v1 = TreePattern(pattern2_v1, quoted_node_names=False)  # create the TreePattern Instance
    solution = pattern2_v1.find_match(tree, None)
    print("solution 2 version 1 is", list(solution))

    pattern2_v2 = """ (sample_1, sample_2) ; """  # comma is used separate sibling nodes
    pattern2_v2 = TreePattern(pattern2_v2, quoted_node_names=False)  # create the TreePattern Instance
    solution = pattern2_v2.find_match(tree, None)
    print("solution 2 version 2 is", list(solution))

    pattern2_v3 = """ ("sample_1", "sample_2") ; """  # comma is used separate sibling nodes
    pattern2_v3 = TreePattern(pattern2_v3, quoted_node_names=True)  # create the TreePattern Instance
    solution = pattern2_v3.find_match(tree, None)
    print("solution 2 version 3 is", list(solution))

    pattern2_v4 = """ ('sample_1', 'sample_2') ; """  # comma is used separate sibling nodes
    pattern2_v4 = TreePattern(pattern2_v4, quoted_node_names=True)  # create the TreePattern Instance
    solution = pattern2_v4.find_match(tree, None)
    print("solution 2 version 4 is", list(solution))

    # List all pattern 2 matches
    solution = list(pattern2_v2.find_match(tree, None, maxhits=None))
    print("All solutions for pattern 2 are:")
    print(solution)

    # Find the total number of pattern2 matches
    solution = len(list(pattern2_v2.find_match(tree, None, maxhits=None)))
    print("The number of solutions for pattern 2 are:")
    print(solution)


if __name__ == "__main__":
    #test()
    test_tutorial()
