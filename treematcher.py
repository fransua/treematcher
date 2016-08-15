#system modules
import re
from itertools import permutations

# third party modules
import ast
import six
from ete3 import PhyloTree, Tree, NCBITaxa

# internal modules
# ...


class TreePatternCache(object):
    def __init__(self, tree):
        """ Creates a cache for attributes that require multiple tree
        traversal when using complex TreePattern queries.

        :param tree: a regular ETE tree instance
         """
        # Initialize cache (add more stuff as needed)
        self.leaves_cache = tree.get_cached_content()
        self.all_node_cache = tree.get_cached_content(leaves_only=False)

    def get_cached_attr(self, attr_name, node, leaves_only=False):
        """
        easy access to cached attributes on trees.

        :param attr_name: any attribute cached in tree nodes (e.g., species,
         name, dist, support, etc.)

        :param node: The pattern tree containing the cache

        :leaves_only: If True, return cached values only from leaves

        :return: cached values for the requested attribute (e.g., Homo sapiens,
         Human, 1.0, etc.)

        """
        #print("USING CACHE")
        cache = self.leaves_cache if leaves_only else self.all_node_cache
        values = [getattr(n, attr_name, None) for n in cache[node]]
        return values

    def get_leaves(self, node):
        return self.leaves_cache[node]

    def get_descendants(self, node):
        return self.all_node_cache[node]


class _FakeCache(object):
    """TreePattern cache emulator."""
    def __init__(self):
        pass

    def get_cached_attr(self, attr_name, node, leaves_only=False):
        """ Helper function to mimic the behaviour of a cache, so functions can
        refer to a cache even when one has not been created, thus simplifying code
        writing. """
        if leaves_only:
            iter_nodes = node.iter_leaves
        else:
            iter_nodes = node.traverse

        values = [getattr(n, attr_name, None) for n in iter_nodes()]
        return values

    def get_leaves(self, node):
        return node.get_leaves()

    def get_descendants(self, node):
        return node.get_descendants()


class PatternSyntax(object):
    def __init__(self):
        # Creates a fake cache to ensure all functions below are functioning
        # event if no real cache is provided
        self.__fake_cache = _FakeCache()
        self.__cache = None

    def __get_cache(self):
        if self.__cache:
            return self.__cache
        else:
            return self.__fake_cache

    def __set_cache(self, value):
        self.__cache = value

    cache = property(__get_cache, __set_cache)

    def leaves(self, target_node):
        return sorted([name for name in self.cache.get_cached_attr(
            'name', target_node, leaves_only=True)])

    def descendants(self, target_node):
        return sorted([name for name in self.cache.get_cached_attr(
            'name', target_node)])

    def species(self, target_node):
        return set([name for name in self.cache.get_cached_attr(
            'species', target_node, leaves_only=True)])

    def contains_species(self, target_node, species_names):
        """
        Shortcut function to find the species at a node and any of it's descendants.
        """
        if isinstance(species_names, six.string_types):
            species_names = set([species_names])
        else:
            species_names = set(species_names)

        found = 0
        for sp in self.cache.get_cached_attr('species', target_node, leaves_only=True):
            if sp in species_names:
                found += 1
        return found == len(species_names)

    def contains_leaves(self, target_node, node_names):
        """ Shortcut function to find if a node contains at least one of the
        node names provided. """

        if isinstance(node_names, six.string_types):
            node_names = set([node_names])
        else:
            node_names = set(node_names)

        found = 0
        for name in self.cache.get_cached_attr('name', target_node, leaves_only=True):
            if name in node_names:
                found += 1
        return found == len(node_names)

    def n_species(self, target_node):
        """ Shortcut function to find the number of species within a node and
        any of it's descendants. """

        species = self.cache.get_cached_attr('species', target_node, leaves_only=True)
        return len(set(species))

    def n_leaves(self, target_node):
        """ Shortcut function to find the number of leaves within a node and any
                of it's descendants. """
        return len(self.cache.get_leaves(target_node))

    def n_duplications(self, target_node):
        """
            Shortcut function to find the number of duplication events at or below a node.
            :param target_node: Node to be evaluated, given as @.
            :return: True if node is a duplication, otherwise False.
        """
        events = self.cache.get_cached_attr('evoltype', target_node)
        return(events.count('D'))

    def n_speciations(self, target_node):
        """
            Shortcut function to find the number of speciation events at or below a node.
        """
        events = self.cache.get_cached_attr('evoltype', target_node)
        return(events.count('S'))


    # TODO: review next function

    def smart_lineage(self, constraint):
        """ Get names instead of tax ids if a string is given before the "in
        @.linage" in a query. Otherwise, returns Taxonomy ids. Function also
        works for constraint that contains something besides the given target
        node (e.g., @.children[0].lineage).

        :param constraint: Internal use.

        :return:  Returns list of lineage tax ids if taxid is searched,
        otherwise returns names in lineage. """

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


class TreePattern(Tree):
    def __str__(self):
        return self.get_ascii(show_internal=True, attributes=["name"])

    def __init__(self, newick=None, format=1, dist=None, support=None,
                 name=None, quoted_node_names=True, syntax=None):
        """ Creates a tree pattern instance that can be used to search within
        other trees.

        :param newick: Path to the file containing the tree or, alternatively,
            the text string containing the same information.
        :param format: subnewick format that is a number between 0 (flexible)
            and 9 (leaf name only), or 100 (topology only).
        :param quoted_node_names: Set to True if node names are quoted with
            stings, otherwise False.
        :syntax: Syntax controller class containing functions to use as
            constraints within the pattern.
        """
        # Load the pattern string as a normal ETE tree, where node names are
        # python expressions
        super(TreePattern, self).__init__(newick, format, dist, support, name, quoted_node_names)

        # Set a default syntax controller if a custom one is not provided
        self.syntax = syntax if syntax else PatternSyntax()

    # FUNCTIONS EXPOSED TO USERS START HERE
    def match(self, node, cache=None):
        """
        Check all constraints interactively on the target node.

        :param node: A tree (node) to be searched for a given pattern.

        :param local_vars:  Dictionary of TreePattern class variables and
        functions for constraint evaluation.

        :return: True if a match has been found, otherwise False.
        """
        self.syntax.cache = cache
        # does the target node match the root node of the pattern?
        status = self.is_local_match(node, cache)

        # if so, continues evaluating children pattern nodes against target node
        # children
        if status and self.children:

                # print("Now checking Children")
                if len(node.children) >= len(self.children):
                    # If pattern node expects children nodes, tries to find a
                    # combination of target node children that match the pattern
                    for candidate in permutations(node.children):
                        sub_status = True
                        for i in range(len(self.children)):
                            st = self.children[i].match(candidate[i], cache)
                            if st is not None:
                                sub_status &= st
                        status = sub_status
                        if status:
                            break
                else:
                    status = False

        return status

    def relaxed_match(self, node, cache=None):
        """
        Check all constraints interactively on the target node.

        :param node: A tree (node) to be searched for a given pattern.

        :param local_vars:  Dictionary of TreePattern class variables and
        functions for constraint evaluation.

        :return: True if a match has been found, otherwise False.
        """
        self.syntax.cache = cache
        # does the target node match the root node of the pattern?
        status = self.is_local_match(node, cache)
        print("Matching target node {} to pattern node {}?".format(node.name, self.name))
        print("status is {}".format(status))

        # if so, continues evaluating children pattern nodes against target node
        # children
        matched_nodes = []
        if status:
            if self.name != '*':
                matched_nodes += [node]
            if self.children:
                print("node children are {}".format(node.children))

                print("Now checking Children")
                if len(node.children) >= len(self.children):
                    # If pattern node expects children nodes, tries to find a
                    # combination of target node children that match the pattern
                    for candidate in permutations(node.children):
                        sub_status = True
                        sub_nodes = []
                        for i in range(len(self.children)):
                            st, self, sub_sub_nodes = self.children[i].relaxed_match(candidate[i], cache)
                            print("Child status is {}".format(st))
                            if st is not None:
                                sub_status &= st
                            sub_nodes += sub_sub_nodes

                        status = sub_status
                        matched_nodes += sub_nodes
                        if status:
                            break
                else:
                    print("Pattern has more children than target")
                    if self.name == '*':
                        sub_status = True
                        sub_nodes = []
                        for i in range(len(self.children)):
                            print("Checking last pattern node")
                            st, self, sub_sub_nodes = self.children[i].relaxed_match(node, cache)
                            print("Second child status is {}".format(st))
                            sub_status &= st
                            sub_nodes += sub_sub_nodes
                        status = sub_status
                        matched_nodes += sub_nodes
                    else:
                        status = False

        self.syntax.cache = None
        return status, self, matched_nodes

    def is_local_match(self, target_node, cache):
        """ Evaluate if this pattern constraint node matches a target tree node.

        TODO: args doc here...
        """

        # Creates a local scope containing function names, variables and other
        # stuff referred within the pattern expressions. We use Syntax() as a
        # container of those custom functions and shortcuts.

        constraint_scope = {attr_name: getattr(self.syntax, attr_name)
                            for attr_name in dir(self.syntax)}
        constraint_scope.update({"__target_node": target_node})

        if not self.name:
            # empty pattern node should match any target node
            constraint = ''
        elif '@' in self.name:
            # converts references to node itself
            constraint = self.name.replace('@', '__target_node')
        elif '*' == self.name:
            # star pattern node should match any target node
            constraint = ''

        else:
            # if no references to itself, let's assume we search an exact name
            # match (allows using regular newick string as patterns)
            constraint = '__target_node.name == "%s"' % self.name

        try:
            st = eval(constraint, constraint_scope) if constraint else True
            # if isinstance(st, six.string_types):
            #     raise ValueError
            # else:
            #     st = bool(st)  # all sets return true

        except ValueError:
            #raise ValueError("not a boolean result: %s. Check quoted_node_names." %
            #                 (st))
            raise ValueError("not a boolean result: . Check quoted_node_names.")

        except (AttributeError, IndexError) as err:
            raise ValueError('Constraint evaluation failed at %s: %s' %
                             (target_node, err))
        except NameError:
            try:
                # temporary fix. Can not access custom syntax on all nodes. Get it from the root node.
                root_syntax = self.get_tree_root().syntax
                # print("root_syntax is ", root_syntax)
                constraint_scope = {attr_name: getattr(root_syntax, attr_name)
                                    for attr_name in dir(root_syntax)}
                constraint_scope.update({"__target_node": target_node})

                st = eval(constraint, constraint_scope) if constraint else True

                # if isinstance(st, six.string_types):
                #    raise ValueError
                # else:
                #    st = bool(st)  # all sets return true
                return st
            except NameError as err:
                raise NameError('Constraint evaluation failed at %s: %s' %
                         (target_node, err))
        else:
            # self.syntax.__cache = None
            return st

    def find_match(self, tree, maxhits=1, cache=None, target_traversal="preorder", relaxed_match=False):
        """ A pattern search continues until the number of specified matches are
        found.
        :param tree: tree to be searched for a matching pattern.
        :param maxhits: Number of matches to be searched for.
        :param cache: When a cache is provided, preloads target tree information
                               (recommended only for large trees)
        :param None maxhits: Pattern search will continue until all matches are found.
        :param target_traversal: choose the traversal order (e.g. preorder, levelorder, etc)
        :param nested: for each match returned,

        """
        if relaxed_match == False:

            num_hits = 0
            for node in tree.traverse(target_traversal):
                if self.match(node, cache):
                    num_hits += 1
                    yield node
                if maxhits is not None and num_hits >= maxhits:
                    break

        else:
            num_hits = 0
            matched_node_list = []
            for node in tree.traverse(target_traversal):
                status, pattern_node, matched_nodes = self.relaxed_match(node, cache)
                matched_node_list += matched_nodes
                if self != pattern_node:
                    self = pattern_node
                if status:
                    num_hits += 1
                    yield matched_node_list
                if maxhits is not None and num_hits >= maxhits:
                    break



def test1():
    pattern = """ ((((d)c)b)a);"""
    tp = TreePattern(pattern, format=1, quoted_node_names=False)
    print tp


def test():
    '''
    181 leaves tree test
    '''
    import time

    t = PhyloTree(
        "((((((((((((((((Human_1, Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), ((((((((((((((((Human_1, Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1))), ((((((((((((((((Human_1, Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1)), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), ((((((((((((((((Human_1, Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), ((((((((((((((((Human_1, Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1)), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1));", format=8)

    #print len(t)

    t.set_species_naming_function(lambda n: n.name.split("_")[0] if "_" in n.name else '')
    t.get_descendant_evol_events()
    cache = TreePatternCache(t)
    pattern = TreePattern(
        """  (('n_duplications(@) > 0')'n_duplications(@) > 0 ')'contains_species(@, ["Chimp", "Human"])' ; """)

    #basic usage
    start_time = time.time()
    for i in range(0, 1000):
        list(pattern.find_match(t, maxhits=None))
    end_time = time.time()
    total_time= (end_time - start_time) / 1000.00
    print("time without cache", total_time)


    # Using Cache
    start_time_cache = time.time()
    for i in range(0, 1000):
        list(pattern.find_match(t, maxhits=None, cache=cache))
    end_time_cache = time.time()
    total_time_cache = (end_time_cache - start_time_cache) / 1000.00
    print("time with cache", total_time_cache)


def relaxed_example():

    """
    The code currently works for trees 1 and 2, but returns last node in match... want to return all nodes that match?

    Note the pattern is not greedy (matches minimum pattern) and t1 only returns 1 match though there are technically 2

    Now works for 3 and 4 but these return first node of match



    """

    t1 = """ (((((( e ) e ) d ) c ) b ) a) ; """
    t1 = PhyloTree(t1, format=8, quoted_node_names=False)

    t2 = """     ((((((((e) e ) d ) c ) a ) e ) c ) a) ; """
    t2 = PhyloTree(t2, format=8, quoted_node_names=False)

    t3 = """ ( e ) a ; """
    t3 = PhyloTree(t3, format=8, quoted_node_names=False)

    #pattern = TreePattern(""" ((e) * ) a ;""", quoted_node_names=False)
    pattern = TreePattern(""" ((((e ) * ) c ) * ) a ;""", quoted_node_names=False)

    border = "\n" + "#" * 110 + "\n"
    print(border)
    print("searching tree 1")
    print list(pattern.find_match(t1, maxhits=1, relaxed_match=True))

    print(border)
    print("searching tree 2")
    print list(pattern.find_match(t2, maxhits=1,relaxed_match=True))
    print(border)

    print(border)
    print("searching tree 3")
    print list(pattern.find_match(t3, maxhits=1,relaxed_match=True))
    print(border)





if __name__ == "__main__":
    #test()
    relaxed_example()


