# TreeMatcher: A new tool for creating Python-based queries on trees

### Program Description

In mathematics, a standard way of representing graphical trees with edge lengths is the Newick format. The TreeMatcher module extends the Newick format to define a tree pattern and includes rules and filters with a Python-based vocabulary. These patterns are then searched for using a tree traversal algorithm. A pattern can be written by accessing the attributes and functions available to an ETE Tree (see Tree and PhyloTree classes), using Python code directly, or through custom functions and syntax.

### How to use treematcher

The simplest way to begin using treematcher is to create a pattern on a single node. In the following example, a string defines the pattern and a TreePattern instance is created. If an attribute is not specified, the node name is assumed by default.

### Example 1
```
# Find a node named "sample_1"
pattern1 = ' sample_1 ; '	 # begin with a string
pattern1 = TreePattern(pattern1)  # create a TreePattern Instance

```
Now that you know how to search for the name of a single node, you may be tempted to access other nodes through constraints like @.children[0].name=="sample_1" and @.children[1].name=="sample_2" but calling a node's descendants in this way restricts the order in which they are considered a match. For example, the permutation @.children[0].name=="sample_2" and @.children[1].name=="sample_1" would not be returned as a match. Using the Newick format ensures that both permutations of children are matched.

### Example 2
```
# Find a tree where sample_1 and sample_2 are siblings.
pattern2 = TreePattern(' (sample_1, sample_2) ; ')
```

Note that the format type is set to 1 as the default which does not allow internal node names. Access other Newick format types using the format argument.

### Example 3
```
# Find a tree where sample_1 and sample_2 are children of the parent ancestor_a.
pattern3 = ' (sample_1, sample_2) ancestor_a ; '
pattern3 = TreePattern(pattern3, format=8)
```

### Quoted node names and the node symbol @
In order to differentiate the parentheses of a function call from the parentheses defining Newick structure, quoted node names are used. The quotes surrounding each node will be removed and the contents inside will be processed as python code. If no quotes are present, as in the previous examples, all parenthesis are assumed to be part of the Newick structure. In order to ensure that the pattern is being processed correctly, set the quoted_node_names to False when not quoting node names. You can set quoted_node_names to True when they are used but this is assumed by default. In order to access a method on a node, use the @ symbol to represent the node.

```
# pattern tree where sample_1 and sample_2 are siblings.
pattern1 = TreePattern(' (sample_1, sample_2) ; ', quoted_node_names=False)
# pattern tree where sample_2 and another leaf are siblings.
pattern2 = TreePattern(""" ('len(@.children)==0', 'sample_2') ; """, quoted_node_names=True)
# pattern tree where sample_1 and another leaf are siblings.
pattern3 = TreePattern(""" ('sample_1', '@.is_leaf()') ; """, quoted_node_names=True)
```

### To Run
To run, use the find_match function. By default, find_match will look for every match. If one match needs to be searched, set maxhits to 1. To know the number of matches, use len().

```
tree = Tree("((sample_1,sample_2)ancestor_a,(sample_1,sample_2)ancestor_b)root;", format = 8)
# find the parent node of the siblings sample_1 and sample_2
pattern = TreePattern(' (sample_1, sample_2) ; ', quoted_node_names=False)
solution = list(pattern.find_match(tree, None))
print("The number of solutions are: ", len(solution))

```

The following tutorial shows the previous examples in more detail.

### Tutorial 1: Introduction to patterns using treematcher.

```
    tree = Tree("((sample_1,sample_2)ancestor_a,(sample_1,sample_2)ancestor_b)root;", format=8)
    print tree

    border = "\n" + "#" * 110 + "\n"

    #######################################################
    print(border)
    print("Find all nodes named sample_1.")
    print(border)
    #######################################################

    # search for the name attribute.
    # The name is quoted but the node is not so quoted_node_names is set to False.
    pattern1_v1 = """ @.name=="sample_1" ; """
    pattern1_v1 = TreePattern(pattern1_v1, quoted_node_names=False)  # TP Instance
    solution = pattern1_v1.find_match(tree, None)
    print("version 1", list(solution))

    # When no attribute is given, the node name is assumed
    pattern1_v2 = """ sample_1 ; """
    pattern1_v2 = TreePattern(pattern1_v2, quoted_node_names=False)
    solution = pattern1_v2.find_match(tree, None)
    print("version 2", list(solution))

    # Find the total number of pattern2 matches
    solution = len(list(pattern1_v2.find_match(tree, None)))
    print("The number of solutions for pattern 1 is:", solution)

    #######################################################
    print(border)
    print("Find a tree where sample_1 and sample_2 are siblings.")
    print(border)
    #######################################################

    # Only need to find if there is a single match so use maxhits=1
    pattern2_v1 = """ (sample_1, sample_2) ; """  # comma is used separate sibling nodes
    pattern2_v1 = TreePattern(pattern2_v1, quoted_node_names=False)  # create the TreePattern Instance
    solution = pattern2_v1.find_match(tree, maxhits=1)
    print("solution ", list(solution))

    #If you want to know if the match exists at a specific node, use match()
    solution1 = pattern2_v1.match(tree.children[0])
    solution2 = pattern2_v1.match(tree.children[1])
    print("solution 1 is", solution1)
    print("solution 2 is", solution2)

    # If you want all matches, use maxhits value of None.
    all_solutions = list(pattern2_v1.find_match(tree, None))
    print("All solutions for pattern 2 are:", all_solutions)


```

### The results of the Tutorial 1 are as follows:

```

      /-sample_1
   /-|
  |   \-sample_2
--|
  |   /-sample_1
   \-|
      \-sample_2

##############################################################################################################

Find all nodes named sample_1.

##############################################################################################################

('version 1', [Tree node 'sample_1' (0x10ac1365), Tree node 'sample_1' (0x10ac3315)])
('version 2', [Tree node 'sample_1' (0x10ac1365), Tree node 'sample_1' (0x10ac3315)])
('The number of solutions for pattern 1 is:', 2)

##############################################################################################################

Find a tree where sample_1 and sample_2 are siblings.

##############################################################################################################

('solution ', [Tree node 'ancestor_a' (0x108f09c9)])
('solution 1 is', True)
('solution 2 is', True)
('All solutions for pattern 2 are:', [Tree node 'ancestor_a' (0x108f09c9), Tree node 'ancestor_b' (0x1091024d)])
```

A short list of commonly used constraints is given in the following table.

Table 1: Examples of common constraints.

|  type                     |custom|  syntax example       						    | example meaning       				        |  Comments																        |
| --------------------------|:-:|:-------------------------------------------------:|:---------------------------------------------:|:-----------------------------------------------------------------------------:|
| node                      |   | @	            						            |a  node, default for nodes left blank	        | Use @.attribute to access attribute, function(@) to access function           |
| node name                 |   | node_name, "node_name", or @.name=="node_name"	| equivalent to @.name == "sample1" 	        | Looking for multiple names, use list: @.name in ("sample1","sample2")         |
| distance                  |   | @.dist >= 0.5     					            | branch length no less than 0.5		        | Use any of the following: <, <=, ==, >=, !=								    |
| support                   |   | @.support > 0.9	            		            | Has a support value greater than 0.90	        | 																		        |
| species                   |   | @.species=="Homo sapiens"	    		            | Homo sapiens is species of node		        | See set_species_naming_function()	for details							        |
| scientific name           |   | @.sci_name == Euarchontoglires 		            | scientific name is Euarchontoglires	        | See annotate_ncbi_taxa() function for details 						        |
| rank                      |   | @.rank == subfamily 					            | node is ranked at the subfamily level	        | See annotate_ncbi_taxa() function for details							        |
| taxonomic id              |   | @.taxid == 207598						            | 20758 is the taxid of the node			    | See annotate_ncbi_taxa() function for details							        |
| number of children        |   | len(@.children)						            | binary tree internal node has 2, leaf hss 0   | use quoted node names to differentiate parentheses from Newick Structure	    |
| size of subtree           |   | @.size > 5							            | number of descendants/size of tree	        | custom attribute 														        |
| number of leaves          |   | len(@)                                            |                                               |                                                                               |
| lineage                   | * | 9606 in @.lineage or "Homo sapiens" in @.lineage  | Homo sapiens in @.lineage				        | Find NCBI taxonomy ID or the full scientific name	in a node's lineage	        |
| species in descendant node| * | "Homo sapiens" in @.contains_species	            | find species contained in a leaf		        | custom attribute														        |
| leaf name                 | * | Pan_troglodytes_1 in leaves(@)		            | Pan_troglodytes_1 is descendant leaf name	    | Find the leaf name within a list of leaf names, custom treematcher attribute  |
*custom functions and attributes do not exist outside of the treematcher class or have different functionality than found elsewhere in ETE.



### Advanced Topics

#### optimization and using a cache

Virtually any attribute available in ETE can be searched for on a tree, however, the larger the structure the more complex the pattern is, the more computationally intensive the search will be. Large Newick trees with complex conditional statements calling functions that require several tree traversals is not recommended.
Instead, break complex patterns into smaller searches. If conditional statements are used, try putting the part of the search that you think will be faster first.


```
    import time

    t = PhyloTree(
        "((((((((((((((((Human_1, Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1), Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1);", format=8)
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


```

####  Custom Functions
You can use your own custom functions and syntax in treematcher.  In the following example, a custom function is created in a custom class called MySyntax.

```
# Expanding vocabulary
class MySyntax(PatternSyntax):
	def my_nice_function(self, node):
		return node.species == 'Chimp'

my_syntax = MySyntax()

pattern = """ 'my_nice_function(@)'; """
t_pattern = TreePattern(pattern, syntax=my_syntax)
for match in t_pattern.find_match(t, cache):
	print(list(match))

```

### Command line tool
|  argument       						| meaning       						|
| ----------------------------------	|:-------------------------------------:|
| --pattern								| string of  pattern in Newick format	|
| --pattern-format						| format for patter, default = 1		|
| --trees								| list of trees to search				|
| --tree-format							| format for trees, default = 1			|									|
| --quoted-node-names 					| default = False						|


ete3 treematcher --pattern "sample_1, sample_2;" --pattern-format 8 --tree-format 8 --trees "sample_3,(sample_1,sample_2)sample_0;" --quoted-node-names

