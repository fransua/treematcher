import re
import itertools
from collections import defaultdict, deque

import six
import copy
from ete3 import PhyloTree, Tree, NCBITaxa

from symbols import SYMBOL, SET
from pprint import pprint

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

class TreePattern(Tree):
    def __str__(self):
        return self.get_ascii(show_internal=True, attributes=["name"])

    def __repr__(self):
        return "Pattern node '%s' (%s)" %(self.name, hex(self.__hash__()))

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
        self.controller = None

    def is_in_bounds(self, bound, test):
        """
        Checks in case of indirect connection the current node is inside the bounds.
        :param bound: specifies what bound to test
        """

        if bound == 'low':
            return test >= self.controller["low"]
        elif bound == 'high':
            if self.controller["high"] >= 0:
                # high bound has to be less but not equal, because it performs a 'can skip one more' test.
                return test <= self.controller["high"]
        return True

    def find_extreme_case(self, match_list):
        """
        Given a single match pattern and a list of (matched) nodes
        returns the one that fulfills the constraint. That constraint should
        describe an axtreme condition such as minimum or maxhimum.

        :param match_list: a list of (matched) nodes
        """

        # find the node in self that has the comparison expression.
        # now is useless. In case this comparison can be made inside other
        # patterns, the every pattern should compare their corresponding node of
        # the other patterns.
        for node in self.traverse():
            if node.controller["single_match"]:
                pattern = node
                constraint = pattern.controller["single_match_contstraint"]

        # for all nodes in the matched pattern find the one that fits the best
        # the constraint expression.
        # assumes that the root has to be tested.
        correct_node = match_list[0]

        st = False # truth flag for eval()

        for node in match_list:
            # for every node in the list compare the best match so far with
            # the current node in the list priority.
            constraint = pattern.controller["single_match_contstraint"]
            constraint_scope = {attr_name: getattr(self.syntax, attr_name)
                                for attr_name in dir(self.syntax)}
            constraint_scope.update({"__target_node": node})
            constraint_scope.update({"__correct_node": correct_node})
            constraint = constraint.replace("@", "__target_node")
            constraint = constraint.replace(SET["all_nodes"], "__correct_node")

            try:
                st = eval(constraint, constraint_scope)

            except ValueError:
                raise ValueError("not a boolean result: . Check quoted_node_names.")

            except (AttributeError, IndexError) as err:
                raise ValueError('Constraint evaluation failed at %s: %s' %
                                 (target_node, err))
            except NameError:
                try:
                    # temporary fix. Can not access custom syntax on all nodes. Get it from the root node.
                    root_syntax = self.get_tree_root().syntax
                    constraint_scope = {attr_name: getattr(root_syntax, attr_name)
                                        for attr_name in dir(root_syntax)}
                    constraint_scope.update({"__target_node": node})
                    constraint_scope.update({"__correct_node": correct_node})

                    st = eval(constraint, constraint_scope)
                    if st: correct_node = node

                except NameError as err:
                    raise NameError('Constraint evaluation failed at %s: %s' %
                             (target_node, err))
            else:
                if st: correct_node = node
        return correct_node

    def decode_repeat_symbol(self, bounds):
        """
        Extracts valuable information from the controlled skipping case.

        :param bounds: a string that contains the {x-y} pattern.

        returns a list with the lower and the highest bounds in respectively order.
        """
        bounds = bounds[1:len(bounds)-1]
        if '-' in bounds:
            split = bounds.split("-")
            low = int(split[0]) if split[0] else 0
            high = int(split[1] ) if split[1] else -1
        else:
            low = high = int(bounds)
        return [low, high]

    def parse_metacharacters(self, raw_constraint):
        """
        Takes a string as node name and extracts the metacharacters.
        Assumes that all metacharacters are defined at the end of the string.
        returns a list containing all metacharacters in the string.
        """
        metacharacters = []
        while len(raw_constraint) > 0 and  raw_constraint[len(raw_constraint)-1] in SYMBOL.values():

            if SYMBOL["defined_number_set_start"] in raw_constraint:
                metacharacters += [ raw_constraint[
                    raw_constraint.find(SYMBOL["defined_number_set_start"]):raw_constraint.find(SYMBOL["defined_number_set_end"])
                    + 1 ]]
                raw_constraint = raw_constraint[0:raw_constraint.find(SYMBOL["defined_number_set_start"])]
            else:
                metacharacters += [ raw_constraint[len(raw_constraint) - 1] ]

                if len(raw_constraint) > 0:
                    raw_constraint = raw_constraint[0:len(raw_constraint)-1]
        return raw_constraint, metacharacters

    def parse_node_name(self):
        """transforms node.name into an evaluable python expression. Metachars are also
        extracted and parsed.

        """
        clean_name, metachar = self.parse_metacharacters(self.name)

        # Parse metachars
        # ...

        # Translate alias and shortcut expressions in clean names
        if '@' not in clean_name:
            constraint = '__target_node.name == "%s"' %clean_name
        else:
            constraint = self.constraint.replace('@', '__target_node')

        if self.children:
            constraint = '(%s) and __target_node.children' %constraint
        else:
            constraint = '(%s) and not __target_node.children' %constraint
        
        return constraint, metachar

    def define_skipping_properties(self, metacharacters):
        controller = {}

        for metacharacter in metacharacters:
            # update controller in case of root or leaf metacharacters and set the node name
            if metacharacter == SYMBOL["is_root"]:
                controller["root"] = True
            if metacharacter == SYMBOL["is_leaf"]:
                controller["leaf"] = True

            # update controller according to metacharacter connection properties.
            if metacharacter == SYMBOL["one_or_more"]:
                controller["allow_indirect_connection"] = True
                controller["direct_connection_first"] = False
                controller["low"] = 1
            elif metacharacter == SYMBOL["zero_or_more"]:
                controller["allow_indirect_connection"] = True
                controller["direct_connection_first"] = True
            elif metacharacter == SYMBOL["zero_or_one"]:
                controller["direct_connection_first"] = True
                controller["allow_indirect_connection"] = True
                controller["high"] = 1
            elif SYMBOL["defined_number_set_start"] in metacharacter:
                split = self.name.split(SYMBOL["defined_number_set_start"])
                bounds = self.decode_repeat_symbol(metacharacter)
                controller["low"] = bounds[0]
                controller["high"] = bounds[1]
                controller["allow_indirect_connection"] = True
                if controller["low"]  == 0: controller["direct_connection_first"] = True
                else: controller["direct_connection_first"] = False

        return controller

    def init_controller(self):
        """
        Creates a dictionary that contains information about a node.
        That information is about how a node interacts with the tree topology.
        It describes how the metacharacter connects with the rest of nodes and
        if it is leaf or root.
        """

        # Interpret node name to python expression
        self.constraint, metachars = self.parse_node_name()
        # bounds
        self.controller = {}
        self.controller["min_instances"] = 0
        self.controller["max_instances"] = -1

        # review scope
        # transform sets to the corresponding code
        #if SET["any_child"] in self.name:
        #    self.name = " any( " + self.name.split("[")[0] + " " + ("[" + self.name.split("[")[1]).replace(SET["any_child"], "x") + " for x in __target_node.children)"
        #if SET["children"] in self.name:
        #    self.name = " all( " + self.name.split("[")[0] + " " + ("[" + self.name.split("[")[1]).replace(SET["children"], "x") + " for x in __target_node.children)"
        #if SET["all_nodes"] in self.name:
        #    controller["single_match"] = True
        #    controller["single_match_contstraint"] = self.name
        #    self.name = '@'

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

        #check the zero intermediate node case.
        #assumes that SYMBOL["zero_or_more"] has only one child.
        if self.controller["direct_connection_first"] and not self.is_leaf():
            self = self.children[0]

        status = self.is_local_match(node, cache)

        if not status:
            if self.controller["allow_indirect_connection"] and self.is_leaf():
                pass
            elif self.up is not None and self.up.controller["allow_indirect_connection"] and self.up.is_in_bounds("high", self.up.controller["skipped"] + 1 ):  # skip node by resetting pattern
                status = True
                self = self.up
                self.controller["skipped"] += 1
        elif self.controller["allow_indirect_connection"] and self.controller["skipped"] == 0:
            self.controller["skipped"] += 1


        # if so, continues evaluating children pattern nodes against target node
        # children
        if status and self.children:

                #if the number of children do not match, find where they do and check that
                nodes = []

                if len(node.children) < len(self.children):
                    if self.controller["allow_indirect_connection"]:
                        count = 0
                        for skip_to_node in node.traverse(strategy="levelorder"):
                            # skip to node with correct number of children
                            if len(skip_to_node.children) >= len(self.children):
                                count += 1
                                nodes += [skip_to_node]
                                sisters = skip_to_node.get_sisters()
                                if len(sisters) > 0:
                                   for sister in sisters:
                                       nodes += [sister]

                                break
                        if count < 1:
                            status = False

                    else:
                        #print("setting status to false")
                        status = False

                # If pattern node expects children nodes, tries to find a
                # combination of target node children that match the pattern

                if len(nodes) == 0:
                    nodes = [node]

                for node in nodes:
                    sub_status_count = 0
                    if len(node.children) >= len(self.children):
                        test_children_properties = False
                        continious_matched = []
                        for candidate in permutations(node.children):
                            sub_status = True
                            current_matched = 0

                            i = 0
                            j = 0
                            while i < len(self.children):
                                st = self.children[i].match(candidate[j], cache)

                                #print "testing: " + self.children[i].name + " -- " + candidate[j].name + " -- st: " + str(st) + " -- subst: " + str(sub_status)

                                if self.children[i].is_leaf() and self.children[i].controller["allow_indirect_connection"] and sub_status: # REVIEW ME
                                    test_children_properties = True
                                    if st:
                                        current_matched += 1
                                    j += 1
                                    if j < len(candidate):
                                        continue
                                    else:
                                        #print "updating " + str(current_matched) + " to " + str(continious_matched)
                                        if current_matched > 0: continious_matched += [current_matched]
                                        current_matched = 0
                                else :
                                    if st and not self.is_in_bounds("low", self.controller["skipped"]):
                                        # in case it matches, but has exited the lower bound (in case it exist), force False the match
                                        st = False
                                    if st == False and self.controller["allow_indirect_connection"] and len(candidate[i].children) > 0:
                                        pass
                                    else:
                                        sub_status &= st
                                        if sub_status: sub_status_count += 1
                                i += 1
                                j += 1

                            if test_children_properties:
                                continue
                            elif sub_status and sub_status_count > 0:
                                status = True
                                break
                            else:
                                status = False
                        if test_children_properties:
                            sub_sub_status = False
                            #print continious_matched
                            if len(continious_matched) > 0 :
                                for con_matches in continious_matched:
                                    low_test = self.children[-1].is_in_bounds("low", con_matches)
                                    high_test = self.children[-1].is_in_bounds("high", con_matches)
                                    sub_sub_status |= low_test and high_test
                                sub_status = sub_sub_status
                            #print "exited with sub_sub_status: " + str(sub_status)
                            #if not self.children[-1].controller["direct_connection_first"]:
                            elif self.children[-1].controller["direct_connection_first"]:
                                sub_sub_status = True
                            sub_status = sub_sub_status
                            status = sub_status
                            if sub_status: sub_status_count += 1

                    if status and sub_status_count > 0:
                        break

        # 'skipped' tracks the maximum skipped node. So only in case of not match, it decreases
        if not status and self.controller["allow_indirect_connection"]: self.controller["skipped"] -= 1
        return status


    def is_local_match(self, target_node, cache):
        """ Evaluate if a tree nodes matches the constraints in this pattern node.  """

        # Creates a local scope containing function names, variables and other
        # stuff referred within the pattern expressions. We use Syntax() as a
        # container of those custom functions and shortcuts.
        constraint_scope = {attr_name: getattr(self.syntax, attr_name)
                            for attr_name in dir(self.syntax)}
        constraint_scope.update({"__target_node": target_node})

        try:
            if self.constraint:
                st = eval(self.constraint, constraint_scope)
            else:
                st = True

        except ValueError:
            raise ValueError("not a boolean result: . Check quoted_node_names.")

        except (AttributeError, IndexError) as err:
            raise ValueError('Constraint evaluation failed at %s: %s' %
                             (target_node, err))
        except NameError:
            try:
                # temporary fix. Can not access custom syntax on all nodes. Get it from the root node.
                root_syntax = self.get_tree_root().syntax
                constraint_scope = {attr_name: getattr(root_syntax, attr_name)
                                    for attr_name in dir(root_syntax)}
                constraint_scope.update({"__target_node": target_node})

                st = True
                for constraint in constraints:
                    if constraint:
                        st &= eval(constraint, constraint_scope)
                    else: st &= True
                return st
            except NameError as err:
                raise NameError('Constraint evaluation failed at %s: %s' %
                         (target_node, err))
        else:
            return st


    def find_match(self, tree, maxhits=1, cache=None, target_traversal="preorder"):
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

        # the controller is for only one use, because it changes the node names.
        # for that reason a deepcopy of the pattern and not the pattern
        # is traversed. Very slow operation.
        one_use_pattern = copy.deepcopy(self)

        # indicaes whether the pattern is a pattern of single match
        # such as maximum or minimum value of an attribute
        single_match_pattern = False

        # sets the controller for every node.
        # should chagne if only is a metacharacter acording to benchmarking tests.
        # if one node indicates single_match, the whole pattern is single_match.
        # no functionality to support other than pattern's root single_match yet
        for node in one_use_pattern.traverse():
            single_match_pattern |= node.set_controller()

        # in case of single_match pattern save the nodes and filter them
        # in other case find the requested match
        if single_match_pattern:
            matched_nodes = []
            for node in tree.traverse(target_traversal):
                if one_use_pattern.match(node, cache):
                    matched_nodes += [node]
            yield one_use_pattern.find_extreme_case(matched_nodes)
        else:
            num_hits = 0
            for node in tree.traverse(target_traversal):
                if one_use_pattern.match(node, cache):
                    num_hits += 1
                    yield node
                if maxhits is not None and num_hits >= maxhits:
                    break


def compute_match_matrix(pattern, tree):
    c2nodes = defaultdict(set)
    for n in tree.traverse():
        for cn in pattern.traverse():
            if cn.is_local_match(n, None):
                c2nodes[cn.constraint].add(n)
    for k,v in c2nodes.items():
        print "%s\n   %s" %(k, v)
    return c2nodes

def children_match_2(tnode, pnode, c2nodes):
    # If not children expected in pattern node, return True, as local
    # conditions have already been checked
    if not pnode.children:
        yield True

    # We generate all possible combinations of tree_node children and
    # pattern_node children. At least one of those combination should provide all true. 
    for combination in [zip(x, pnode.children)
                        for x in itertools.permutations(tnode.children, len(pnode.children))]:
        match = combination
        for node_ch, pattern_ch in combination:
            if node_ch not in c2nodes[pattern_ch.constraint]:
                match = None
                break
            else:
                if not children_match(node_ch, pattern_ch, c2nodes):
                    match = None
                    break
        if match:
            yield True

def children_match(tnode, pnode, c2nodes):
    # If not children expected in pattern node, return True, as local
    # conditions have already been checked
    if not pnode.children:
        return True

    t_children = set(tnode.children)

    # CHANGE THIS
    min_occur = 0
    max_occur = 2

    matches = []
    matched_children = set()
    for pnode_ch in pnode.children:
        match_nodes = c2nodes[pnode_ch.constraint] & t_children
        if len(match_nodes) > max_occur:
            return False

        # at least a node matching each pattern constraint
        if not match_nodes and min_occur > 0:
            return False

        # Record all children nodes with matches
        matched_children.update(match_nodes)

        # And prepare all permutations of matching in node with minimum occurrences
        if not match_nodes and min_occur == 0:
            matches.append([[None]])
        else:
            matches.append(list(itertools.permutations(match_nodes, max(min_occur, 1))))

    # there should be nodes without a match
    if len(matched_children) < len(t_children):
        return False

    # Let's check if there is a non-overlapping combination of nodes matches
    # satisfies patterns. For instance, avoid cases where one node matches the
    # two required patterns
    for comb in itertools.product(*matches):
        valid = set()
        match = comb
        for x in comb:
            x = set(x)
            inter =  valid & set(x)
            if not inter:
                valid.update(x)
            else:
                match = None
                # let's check next comb
                break

        if match:
            # Let's check inside the node
            for i, pnode_ch in enumerate(pnode.children):

                for node_ch in match[i]:
                    if node_ch is None:
                        continue
                    if not children_match(node_ch, pnode_ch, c2nodes):
                        return False 
            return True

    return False

def expand_loose_connection_aliases(nw):
    def find_first_unmatched_closing_par(string):
        open_par = 0
        for i, ch in enumerate(string):
            if ch == '(':
                open_par += 1
            elif ch == ')':
                if open_par == 0:
                    return i
                open_par -= 1
        return -1
    chunks = []
    i = 0
    while True:
        next_split = nw[i:].find('^')
        if next_split == -1:
            chunks.append(nw[i:])
            break

        chunks.append(nw[i:next_split])
        chunks.append(',')

        i += next_split + 1
        ipoint = find_first_unmatched_closing_par(nw[i:])
        chunks.append(nw[i:i+ipoint+1])
        chunks.append('^')
        i += ipoint + 1

    return ''.join(chunks)

def loose_connection_nodes(tree):
    nodes = []
    for n in tree.traverse():
        if n.name.startswith('^'):
            nodes.append(n)
    return nodes


def find_matches(tree, pattern):
    #for n in split_by_loose_connection(pattern):
    #    print n

    for n in pattern.traverse():
        n.init_controller()

    c2nodes = compute_match_matrix(pattern, tree)
    for root in c2nodes[p1.constraint]:
        if children_match(root, p1, c2nodes):
            print "MATCH"
            print root
    

def test():
    #nw = expand_loose_connection_aliases('(  ( (A,B) ^ (C,D) ) ^ C);')

    t1 = Tree("(((A, (B,C,D)), ((B,C), A)), F);")
    p1 = TreePattern("((C,B,D), A);")



    
    #t1 = Tree("(a, b, b, a);")
    #p1 = TreePattern("(a, a, b);")


    #t1 = Tree("((a, b), c);")
    #p1 = TreePattern("((a, b, d), c);")

    #t1 = Tree("(  (((B,H), (B,B,H), C), A), (K, J));")
    #p1 = TreePattern("((C, (B,H), (B,H)), A);")

    find_matches(t1, p1)    
if __name__ == '__main__':
    test()
