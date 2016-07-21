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
        ("__target.leaves",  "get_cached_attr('name', __target)"),
        ("__target.size", "len(get_cached_attr('name', __target))"),
        ("__target.contains_species", "get_cached_attr('species', __target)"),

        #("__target.size", "len(temp_leaf_cache[__target])"),
        # ("__target.contains_species", "[n.species for n in temp_leaf_cache[__target]]"),
        #("__target.leaves", "[n.name for n in temp_leaf_cache[__target]]"),

    ]

    def __init__(self, *args, **kargs):
        kargs["format"] = 1

        is_exact = False

        self.cache_flag = 0
        self.all_node_cache = None

        self.evol_events_flag = 0
        self.evol_events = None


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

                for search_type in cached_searches:
                    if pattern_string.find(search_type) != -1:
                        self.cache_flag = 1
                if pattern_string.find("is_duplication")!= -1 or pattern_string.find("is_speciation") != -1:
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


    def constrain_match(self, __target, local_vars=None):
        """
        Evaluate constraint on a single node. This functions checks whether a single node matches a single constraint.
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
        Check all constraints on a tree. Permutes the tree and checks that all constraints return True.
        Funtion checks whether a single node and its children match a given pattern.
        :param node: A tree (root node) to be searched for a given pattern.
        :param local_vars:  Dictionary of TreePattern class variables and functions for constraint evaluation.
        :return: True is a match has been found, otherwise False.
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
        Creates cache for attributes that require multiple tree traversal.
         Cache is available on all nodes of the tree.
        :param tree: Pattern being searched for.
        """
        if self.evol_events_flag == 1:
            self.evol_events = tree.get_descendant_evol_events()
            print("self.evol_events are first", self.evol_events)
        #self.temp_leaf_cache = tree.get_cached_content()

        self.all_node_cache = tree.get_cached_content(leaves_only=False)

        # pass by reference, so dictionary is available on all nodes
        for node in self.traverse():
            #node.temp_leaf_cache = self.temp_leaf_cache
            node.all_node_cache = self.all_node_cache
            node.evol_events = self.evol_events

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

            #Search for custom functions
            if self._extra_functions is not None:
                for custom_function in self._extra_functions.keys():
                    node.constraint = re.sub(custom_function + "\(__target",
                                             custom_function + "(__target, pattern", node.constraint)

            #custom function of stuff tpat=whatever

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

    cached_species = __pattern.get_cached_attr('species', __target)
    result = all(sp in cached_species for sp in species_list)
    return result


def contains(__target, __pattern, name_list):
    """
    Shortcut function to find the name(s) within a node or any of it's descendants.
    :param __target: Internal use
    :param __pattern: Internal use
    :param name_list: List of names to be searched for at a node.
    :return: True if all names in name_list are found, otherwise False.
    """
    __pattern.preprocess(__target)

    if isinstance(name_list, six.string_types):
        name_list = [name_list]

    print(pattern0.find_match(tree, None))
    print(pattern1.find_match(tree, None))

def test2():
    tree = PhyloTree(
            "((9315.ENSMEUP00000008285:0.899711,9258.ENSOANP00000027752:0.559777)0.99985:0.11989,((9739.ENSTTRP00000010720:0.164873,9913.ENSBTAP00000003500:0.298158)0.99985:0.109903,((9685.ENSFCAP00000006440:0.239731,(9615.ENSCAFP00000042310:0.122399,(9646.ENSAMEP00000002314:0.18278,9669.ENSMPUP00000005544:0.270727)0.6117:0.0396991)0.99985:0.0702148)0.99985:0.082488,(132908.ENSPVAP00000014833:0.488081,(9796.ENSECAP00000022144:0.310699,(((9785.ENSLAFP00000009512:0.187095,9813.ENSPCAP00000004417:0.493329)0.99985:0.359095,(30611.ENSOGAP00000016876:0.334272,(9483.ENSCJAP00000021314:0.178043,(9601.ENSPPYP00000003401:0.0415077,((61853.ENSNLEP00000003253:0.196659,9544.ENSMMUP00000037769:0.326984)0.835225:0.0989423,(9593.ENSGGOP00000004740:0.101826,9606.ENSP00000182290:0.0204981)0.997196:0.020731)0.307827:0.0046059)0.99985:0.0991112)0.99985:0.162323)0.972253:0.0380139)0.70642:0.0193389,((10141.ENSCPOP00000016274:0.272126,43179.ENSSTOP00000015376:0.458416)0.996119:0.0901785,(37347.ENSTBEP00000013312:0.328061,(10020.ENSDORP00000010739:0.398341,(10116.ENSRNOP00000051746:0.0455948,10090.ENSMUSP00000009396:0.0811741)0.99985:0.269525)0.791467:0.0577236)0.536676:0.0461933)0.99985:0.0620583)0.99985:0.0788824)0.969465:0.0395994)0.635969:0.0171601)0.702925:0.0283261)0.99985:0.11989);",
            format=1, quoted_node_names=False)

    tree.set_species_naming_function(lambda n: n.name.split(".")[0] if "." in n.name else '')
    tree.annotate_ncbi_taxa()

    pattern0 = """
    ( '9443 in @.lineage  and "Primates" in @.lineage and @.name!=9606 ' )' @.support >= 0.9 ';
    """
    pattern0 = TreePattern(pattern0, format=8, quoted_node_names=True)
    #print(len(list(pattern0.find_match(tree, None, maxhits=None))))

    pattern1 = """
    "Mammalia" in @.children[0].lineage  ;
    """
    pattern1 = TreePattern(pattern1, format=8, quoted_node_names=True)
    print(len(list(pattern1.find_match(tree, None, maxhits=3))))

    tree = PhyloTree(
        "((((Anolis_carolinensis_1:1, Gallus_gallus_1:1), (Felis_catus_1:1, (Homo_sapiens_1:1, Pan_troglodytes_1:1)primates)primates), ((Danio_rerio_1:1, (Xenopus_laevis_1:1, Anolis_carolinensis_1:1)), Saccharomyces_cerevisiae_2:1)), Saccharomyces_cerevisiae_1:1);",
        format=1)
    tree.set_species_naming_function(lambda node: node.name.split("_")[0])

    print(len(list(pattern0.find_match(tree, None, maxhits=None))))
    print(len(list(pattern0.find_match(tree, None))))
    print(len(list(pattern1.find_match(tree, None))))
    print(len(list(pattern1.find_match(tree, None, maxhits=None))))

if __name__ == "__main__":
    #test()
    test2()
