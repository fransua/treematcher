# TreeMatcher: A new tool for creating Python-based queries on trees

### Program Description

In mathematics, a standard way of representing graphical trees with edge lengths is the Newick format. The TreeMatcher module defines a tree pattern by extending the Newick format to include rules and filters with a regular-expression-like vocabulary. These patterns are then searched for using a tree traversal algorithm. A pattern can be written by accessing the attributes and functions available to an ETE Tree (see Tree class), using Python code directly, or through custom functions and syntax.

### How to use treematcher

The simplest way to begin using treematcher is to create a pattern to find the name of a single node. In the following example, a string defines the pattern and a TreePattern instance is created. If an attribute is not specified, it is assumed by default that you are searching for a node name.

```
# Find a node named "sample_1"
pattern1 = """ sample_1 ; """	 # begin with a string
pattern1 = TreePattern(pattern1)  # create a TreePattern Instance

```
Now that you know how to search a single node, you may be tempted to access other nodes through constraints like @.children[0].name=="sample_1" and @.children[1].name=="sample_2" but calling a node's descendants in this way restricts the order in which they are considered a match. For example, the equivalent permutation @.children[0].name=="sample_2" and @.children[1].name=="sample_1" would not be returned as a match. Using the Newick format to ensure that both permutations of children are matched.

```
# Find a tree where sample_1 and sample_2 are siblings.
pattern1_v2 = TreePattern(' (sample_1, sample_2) ; ')
```

Note that the format type is set to 1 as the default which does not allow internal node names. Access other Newick format types using the format argument.

```
# Find a tree where sample_1 and sample_2 are children of the parent ancestor_a.
pattern3 = """ (sample_1, sample_2) ancestor_a ; """
pattern3 = TreePattern(pattern3, format=8)
```

### To Run
To run, use the find_match function. By default, find_match will look for every match. If if one match needs to be searched, set maxhits to 1. To know the number of matches, use len().

```
tree = Tree("((sample_1,sample_2)ancestor_a,(sample_1,sample_2)ancestor_b)root;", format = 8)
# find the parent node of the siblings sample_1 and sample_2
pattern = TreePattern(' (sample_1, sample_2) ; ')
solution = list(pattern.find_match(tree, None))
print("The number of solutions are: ", len(solution))

```

See the tutorial at the end of this file for more detailed examples.

##  Examples of common search types


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
| leaf name                 | * | Pan_troglodytes_1 in @.leaves		                | Pan_troglodytes_1 is descendant leaf name	    | Find the leaf name within a list of leaf names, custom treematcher attribute  |
*custom functions and attributes do not exist outside of the treematcher class or have different functionality than found elsewhere in ETE.

To differentiate the parentheses of a set, tuple, or function in your pattern as separate from the parentheses of the Newick format, use quoted node names.


## Advanced Topics

optimization and cacheing

Virtually any attribute available in ETE can be searched for on a tree, however, the larger the structure the more complex the pattern is, the more computationally intensive the search will be. Large Newick trees with complex conditional statements calling functions that require several tree traversals is not recommended.
Instead, break complex patterns into smaller searches. If conditional statements are used, the part of the search that you think will be faster should come first. For example, matching a node name  if faster than finding the distance to all other nodes or finding a name in a lineage. If many tree traversals are required, using a cache is suggested.

```
#cache example
t = PhyloTree(
    "((((Human_1, Chimp_1), (Human_2, (Chimp_2, Chimp_3))), ((Fish_1, (Human_3, Fish_3)), Yeast_2)), Yeast_1);")

t.set_species_naming_function(lambda node: node.name.split("_")[0])
t.get_descendant_evol_events()

#Basic usage
pattern = TreePattern("""(' "Chimp" in species(@)', ''); """)
print pattern.match(t)

#Using cache
cache = TreePatternCache(t)
pattern = TreePattern("""('"Chimp" in species(@)', ''); """)
print pattern.match(t, cache)




```


Custom Functions
Using your own custom functions and syntax is made possible with treematcher. You can alter the function names using a dictionary. For example, suppose you have written two functions: num_of_species and num_of_leaves. Access these functions as follows.


```
pattern4 = """ 'sn(@, 2) and ln(@,2)' ; """

pattern4 = TreePattern(pattern4, format=8, quoted_node_names=True,
                      functions={"sn": number_of_species,
                                 "ln": number_of_leaves})

solution = pattern4.find_match(tree, None, maxhits=None)
print(list(solution))

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
    pattern1_v1 = TreePattern(pattern1_v1, quoted_node_names=False)  # create a TreePattern Instance
    solution = pattern1_v1.find_match(tree, None)
    print("version 1", list(solution))

    # When no attribute is given, the node name is assumed
    pattern1_v2 = """ sample_1 ; """
    pattern1_v2 = TreePattern(pattern1_v2, quoted_node_names=False)  # create a TreePattern Instance
    # print("Pattern 1 version 2 is", pattern1_v2)
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

    # When you only need to find if there is a single match, use maxhits=1
    pattern2_v1 = """ (sample_1, sample_2) ; """  # comma is used separate sibling nodes
    pattern2_v1 = TreePattern(pattern2_v1, quoted_node_names=False)  # create the TreePattern Instance
    solution = pattern2_v1.find_match(tree, maxhits=1)
    print("solution ", list(solution))

    #If you want to know if the match exists at a specific node, use match()
    solution1 = pattern2_v1.match(tree.children[0])
    solution2 = pattern2_v1.match(tree.children[1])
    print("solution 1 is", solution1)
    print("solution 2 is", solution2)

    # If you want all the matches, use maxhits value of None.
    all_solutions = list(pattern2_v1.find_match(tree, None))
    print("All solutions for pattern 2 are:", all_solutions)

```
### The results of the Tutorial 1 are as follows:

`      /-sample_1
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