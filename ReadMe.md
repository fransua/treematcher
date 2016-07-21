# TreeMatcher: A new tool for creating regular-expression-like queries on trees

In mathematics, a standard way of representing graphical trees with edge lengths is the Newick format which uses a specific syntax (such as parentheses and commas). The TreeMatcher module defines a tree pattern by extending the Newick format to include rules and filters with a regular-expression-like vocabulary. These patterns are then searched for using a tree traversal algorithm.

### Syntax

ETE uses its own formalism to represent phylogenetic trees (see Tree class). A pattern can be written by accessing the attributes and functions available to an ETE Tree, using Python code directly, or through custom functions to use as constrinats. A tree pattern structure is created using the Newick format.


The simplest way to write a pattern is to begin with a string. Then use it to create TreePattern instance. In the simplest case, you want to know something about a single node. If an attribute type is not specified, it is assumed that you are searching for a node name by default.

Example 1: Find a node named "sample_1".
```
pattern1 =' sample_1 ;'	 # begin with a string
pattern1 = TreePattern(pattern1)  # create a TreePattern Instance

```
Example 2: Find a tree where sample_1 and sample_2 are siblings.
```
pattern2 = ' sample_1, sample_2 ; '  #  comma is used separate sibling nodes
pattern2 = TreePattern(pattern2)  # create the TreePattern Instance
```

For a short list of commonly used constraints, use the following table.

Table 1: Examples of common constraints.

|  syntax       						| meaning       						|  Comments																|
| ----------------------------------	|:----------------------------------:	|:---------------------------------------------------------------------:|
| @	            						| node, default for nodes left blank	| Use @.attribute to access an attribute, function(@) to access function|
| sample1								| equivalent to @.name == "sample1" 	| Looking for multiple names, use list: @.name in ("sample1","sample2") |
| @.dist >= 0.5     					| branch length no less than 0.5		| 	Use any of the following: <, <=, ==, >=								|
| @.support > 0.9	            		| Has a support value greater than 0.90	| 																		|
| @.species=="Homo sapiens"	    		| Homo sapiens is species of node		| See set_species_naming_function()	for details							|
| 9606 in @.lineage	            		| Homo sapiens in @.lineage				| Find NCBI taxonomy ID or the full scientific name	in a node's lineage	|
| @.sci_name == Euarchontoglires 		| scientific name is Euarchontoglires	| See annotate_ncbi_taxa() function for details 						|
| @.rank == subfamily 					| node is ranked at the subfamily level	| See annotate_ncbi_taxa() function for details							|
| @.taxid == 207598						|20758 is the taxid	of the node			| See annotate_ncbi_taxa() function for details							|
| Pan_troglodytes_1 in @.leaves			| Pan_troglodytes_1 is descendant leaf	| Find the leaf name within a list of leaf names						|
| len(@.children)						| number of children					| use quotes to differentiate parentheses from Newick Structure			|
| "Homo sapiens" in @.contains_species	| find species contained in a leaf		| custom attribute														|
| @.size > 5							| number of descendants/size of tree	| custom attribute 														|
| H_saps_1 in @.leaves					| find name H_saps_1 in leaves			| custom attribute														|


### How to use Newick format to access the structure of a tree
Create a tree structure using Newick format. Replace the name of each node you want to search for with a constraint.


Example 3: Find a tree where sample_1 and sample_2 are children of the parent sample_0.  Note that the format type is set to 1 as the default which does not allow internal node names. Access other Newick format types by setting the format argument.
```
pattern3 = """ (sample_1, sample_2) sample_0 ; """
pattern3 = TreePattern(pattern3, format=8)
```

### To Run
To run, use the find_match function.  By default, find_match will look for one match. If you want to find every match on a tree, set the maximum number of hits to None.


Example 4: For the following tree, find the node that matches pattern2.

tree = Tree("(sample_1,(sample_1,sample_2)sample_0)sample_0;", format = 8)

```
solution = pattern2.find_match(tree, None)
print(list(solution))
```

Example 5: Find the total number of pattern3 matches in the same tree as above.
```
solution = len(list(pattern3.find_match(tree, None, maxhits=None)))
print(solution)
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
