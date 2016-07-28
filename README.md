# TreeMatcher: A new tool for creating regular-expression-like queries on trees

In mathematics, a standard way of representing graphical trees with edge lengths is the Newick format which uses a specific syntax (such as parentheses and commas). The TreeMatcher module defines a tree pattern by extending the Newick format to include rules and filters with a regular-expression-like vocabulary. These patterns are then searched for using a tree traversal algorithm.

### Syntax

ETE uses its own formalism to represent phylogenetic trees (see Tree class). A pattern can be written by accessing the attributes and functions available to an ETE Tree, using Python code directly, or through custom functions to use as constrinats. A tree pattern structure is created using the Newick format.


The simplest way to write a pattern is to begin with a string surrounded by three double quotes at each end. Use the string to create TreePattern instance. In the simplest case, you want to know something about a single node. If an attribute type is not specified, it is assumed that you are searching for a node name by default.

Example 1: Find a node named "sample_1".
```
pattern1 =' sample_1 ;'	 # begin with a string
pattern1 = TreePattern(pattern1)  # create a TreePattern Instance

```
Example 2: Find a tree where sample_1 and sample_2 are siblings.
```
pattern2 = ' (sample_1, sample_2) ; '  #  comma is used separate sibling nodes. Solution is the parent node.
pattern2 = TreePattern(pattern2)  # create the TreePattern Instance
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
| leaf name                 | * | Pan_troglodytes_1 in @.leaves		                | Pan_troglodytes_1 is descendant leaf name	    | Find the leaf name within a list of leaf names, custom treematcher attribute  |
*custom functions and attributes do not exist outside of the treematcher class or have different functionality than found elsewhere in ETE.

### Newick format to access the structure of a tree
Create a tree structure using Newick format. Replace the name of each node you want to search for with a constraint.


Example 3: Find a tree where sample_1 and sample_2 are children of the parent ancestor_a.  Note that the format type is set to 1 as the default which does not allow internal node names. Access other Newick format types by setting the format argument.
```
pattern3 = """ (sample_1, sample_2) ancestor_a ; """
pattern3 = TreePattern(pattern3, format=8)
```

### To Run
To run, use the find_match function.  By default, find_match will look for one match. If you want to find every match on a tree, set the maximum number of hits to None. If you only want to know the number of matches, use len().


Example 4: For the following tree, find the parent node of the siblings sample_1 and sample_2.

tree = Tree("((sample_1,sample_2)ancestor_a,(sample_1,sample_2)ancestor_b)root;", format = 8)

```
solution = list(pattern2.find_match(tree, None))
```

Example 5: Find the total number of patterns that match
```
solution = len(list(pattern2.find_match(tree, None, maxhits=None)))
```


### Custom Functions
Write your own functions and provide them as local variables to the treematcher program. To differentiate the parentheses of a set, tuple, or function in your pattern as separate from the parentheses of the Newick format, use quoted node names. You can alter the function names using a dictionary. For example, suppose you have written two functions: num_of_species and num_of_leaves. Access these functions as follows.

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


Now, let's combine the examples together into a short tutorial.

Tutorial 1: Introduction to patterns using treematcher.

```

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



```
The results of the Tutorial 1 are as follows:

```
      /-sample_1
   /-|
  |   \-sample_2
--|
  |   /-sample_1
   \-|
      \-sample_2

##############################################################################################################

Find a node named sample_1.

##############################################################################################################

(Tree node '@.name=="sample_1"' (0x10b4e081), [Tree node 'sample_1' (0x10b4e0a9)])
(Tree node 'sample_1' (0x10b4e0b5), [Tree node 'sample_1' (0x10b4e0a9)])
(Tree node '"sample_1"' (0x10b4e0ad), [])
(Tree node '"sample_1"' (0x10b4e0bd), [])
All solutions for pattern 1 are:
[Tree node 'sample_1' (0x10b4e0a9), Tree node 'sample_1' (0x10b4e0d5)]
The number of solutions for pattern 1 are:
2
```