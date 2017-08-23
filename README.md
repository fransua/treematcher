# TreeMatcher: A new tool for creating Python-based queries on trees

### Program Description

In mathematics, a standard way of representing graphical trees with edge lengths is the Newick format. The TreeMatcher module extends the Newick format to define a tree pattern and includes rules and filters with a Python-based vocabulary. These patterns are then searched for using a tree traversal algorithm. A pattern can be written by accessing the attributes and functions available to an ETE Tree (see Tree and PhyloTree classes), using Python code directly, or through custom functions and syntax.

### How to use treematcher

The simplest way to begin using treematcher is to create a pattern on a single node. In the following example, a string defines the pattern and a TreePattern instance is created. If an attribute is not specified, the node name is assumed by default.

```
# Example 1: Find a node named "sample_1"
pattern1 = ' sample_1 ; '	 # begin with a string
pattern1 = TreePattern(pattern1)  # create a TreePattern Instance

```
Now that you know how to search for the name of a single node, you may be tempted to access other nodes through constraints like @.children[0].name=="sample_1" and @.children[1].name=="sample_2" but calling a node's descendants in this way restricts the order in which they are considered a match. For example, the permutation @.children[0].name=="sample_2" and @.children[1].name=="sample_1" would not be returned as a match. Using the Newick format ensures that both permutations of children are matched.


Note that the format type is set to 1 as the default which does not allow internal node names. Access other Newick format types using the format argument.

```
# Example 2: Find a tree where sample_1 and sample_2 are siblings under ancestor_a
tree = Tree("((sample_1,sample_2)ancestor_a,(sample_1,sample_2)ancestor_b)root;", format = 8)
pattern2 = TreePattern(' (sample_1, sample_2)ancestor_a ; ', format=8)
```

### Quoted node names and the node symbol @
In order to differentiate the parentheses of a function call from the parentheses defining Newick structure, quoted node names are used. The quotes surrounding each node will be removed and the contents inside will be processed as python code. If no quotes are present, as in the previous examples, all parenthesis are assumed to be part of the Newick structure. In order to ensure that the pattern is being processed correctly, set the quoted_node_names to False when not quoting node names. You can set quoted_node_names to True when they are used but this is assumed by default. In order to access a method on a node, use the @ symbol to represent the node.

```

# Example 3: Find a tree where sample_2 and another leaf are siblings where a leaf is determined by number of children.
pattern2 = TreePattern(""" ('len(@.children)==0', 'sample_2')ancestor_a ; """, quoted_node_names=True)

# Example 4: Find a tree where sample_1 and another leaf are siblings by accessing the the is_leaf method.
pattern3 = TreePattern(""" ('sample_1', '@.is_leaf()')ancestor_a ; """, quoted_node_names=True)
```

### Relax matches

Treematcher allows to test against relax matched patterns.

By setting ` ^ ` as ancestor you enable the loose connection ability.
Loose connection means that either the node(s) under ` ^ ` may connect via any number
of nodes to a upper or higher level node or that the ` ^ ` children may connect
loosely (via any number of nodes) themselves.

```
# example 5: test if tips A, B, C exists in the same tree
pattern5 = TreePattern(" (A, B, C)^ ;")

#example 6: test if (A, B) and (C, D) are conected via any number of nodes
pattern6 = TreePattern(" ((A,B)^), ((C,D)^) ;")

```

You can test for a relax number of children too.
You can use:
` * `, ` + ` and ` { min, max} `
They borrow their meaning from regular expressions.

```
#example 7: zero or more ocuurances of a node
TreePattern( " ('A', '@.dist > 0.5*')'^' ;", quoted_node_names=True)

#example 8: exact number
TreePattern(" ('A', '@.dist > 0,5{2, 5}') ;", quoted_node_names=True)

```


### To Run
To run, use the find_match() function. By default, find_match will look for one match.
To find the number of matches returned, use len().

```

# Example 9: Find the parent node of the siblings sample_1 and sample_2
tree = Tree("((sample_1,sample_2)ancestor_a,(sample_1,sample_2)ancestor_b)root;", format = 8)
pattern = TreePattern(' (sample_1) ; ', quoted_node_names=False)
solution = list(pattern.find_match(tree))
print("The number of solutions are: ", len(solution))

```
The following tutorial shows the previous examples in more detail.

For more details on how to use treematcher read the tutorial.


A short list of commonly used constraints is given in the following table.

Table 1: Examples of common constraints.

|  type                     |custom|  syntax example       						            | example meaning       				        |  Comments																        |
| --------------------------|:-:|:---------------------------------------------------------:|:---------------------------------------------:|:-----------------------------------------------------------------------------:|
| node                      |   | @	            						                    |a  node, default for nodes left blank	        | Use @.attribute to access attribute, function(@) to access function           |
| node name                 |   | node_name, "node_name", or @.name=="node_name"	        | when attribute not specified, name is assumed | Looking for multiple names, use list: @.name in ("sample1","sample2")         |
| distance                  |   | @.dist >= 0.5     					                    | branch length no less than 0.5		        | Use any of the following: <, <=, ==, >=, !=								    |
| support                   |   | @.support > 0.9	            		                    | Has a support value greater than 0.90	        | 																		        |
| species                   |   | @.species=="Homo sapiens"	    		                    | Homo sapiens is species of node		        | See set_species_naming_function()	for details							        |
| scientific name           |   | @.sci_name == Euarchontoglires 		                    | scientific name is Euarchontoglires	        | See annotate_ncbi_taxa() function for details 						        |
| rank                      |   | @.rank == subfamily 					                    | node is ranked at the subfamily level	        | See annotate_ncbi_taxa() function for details							        |
| taxonomic id              |   | @.taxid == 207598						                    | 20758 is the taxid of the node			    | See annotate_ncbi_taxa() function for details							        |
| number of children        |   | len(@.children)						                    | binary tree internal node has 2, leaf hss 0   | use quoted node names to differentiate parentheses from Newick Structure	    |
| size of subtree           |   | @.get_descendants() > 5							        | size of tree  is greater than 5               | number of descendants												            |
| number of leaves          |   | len(@) > 2                                                | number of leaves is greater than 2            | number of leaves descending from a node                                       |
| lineage                   | * | 9606 in @.lineage or "Homo sapiens" in @.named_lineage    | Homo sapiens in @.lineage				        | Find NCBI taxonomy ID or the full scientific name	in a node's lineage	        |
| species in descendant node| * | contains_species(@, ["Pan troglodytes", "Homo sapiens"])	| Find species in the last at or below node     | species at a node and any of it's descendants							        |
| leaf name                 | * | contains_leaves(@, ["Chimp_2", "Chimp_3"])		        | Pan_troglodytes_1 is descendant leaf name	    | Find the leaf name within a list of leaf names                                |
| number of duplications    | * |  		n_duplications(@) > 0                               | Number of duplications beyond and including this node is greater than zero.	    | number of duplication events at or below a node  |
| number of speciations     | * |  		n_speciations(@) > 0                                | Number of speciations beyond and including this node is greater than zero.	    | number of speciation events at or below a node  |

* functions do not exist outside of treematcher classes.



### Advanced Topics

#### optimization

Virtually any attribute available in ETE can be searched for on a tree, however, the larger the structure the more complex the pattern is, the more computationally intensive the search will be. Large Newick trees with complex conditional statements calling functions that require several tree traversals is not recommended.
Instead, break complex patterns into smaller searches. If conditional statements are used, try putting the part of the search that you think will be faster first.

####  Custom Functions
You can use your own custom functions and syntax in treematcher.  In the following example, a custom function is created in a custom class called MySyntax.

```
# Example 8: Expanding vocabulary
class MySyntax(PatternSyntax):
	def my_nice_function(self, node):
		return node.species == 'Chimp'

my_syntax = MySyntax()

pattern = """ 'my_nice_function(@)'; """
t_pattern = TreePattern(pattern, syntax=my_syntax)
for match in t_pattern.find_match(t):
	print(list(match))

```

### Command line tool

|  argument       						| meaning       						                                                  |
| --------------------------------------|:---------------------------------------------------------------------------------------:|
| -p								    | a list of patterns in newick format (filenames with one per file or quoted strings)     |
| -t								    | a list of trees in newick format (filenames or quoted strings)                          |
|-v                                     | prints the current pattern, prints which trees (by number) do not match the pattern     |
| --tree_format							| format for trees, default = 1	                            		                      |
| --quoted_node_names 					| default = True					                            	                      |
| -o, --output                  | output file for search results
| --src_tree_list                       | path to a file containing many target trees, one per line                               |
| --pattern_tree_list                   | path to a file containing many pattern trees, one per line                              |
|-r, --root                             | flag to return the root of the tree if at least a match was found
| --render                              | filename (.SVG, .PDF, or .PNG), to render the tree image                                |
| --tab                                 | output results in tab delimited format, default if -o used and ascii not specified      |
| --ascii                               | output results in ascii format                                                          |

examples:
Read patterns from a file called MyPatterns.txt and apply to each tree in MyTargetTrees.txt, output the results of each pattern in separate files called treematches0.txt, treematches1.txt, etc
If there is only one pattern, the result file will not be numbered.

` ete_search --pattern_tree_list "MyPattern.txt" --tree_format 8 --src_tree_list "MyTargetTrees.txt" -o treematches.txt `

Provide the pattern and tree as strings and print the result to the terminal.
`ete_search -p "(e,d);" --tree_format 8 -t "(c,(d,e)b)a;" `


Count how many trees matches a pattern from a list of trees.
` ete_search -p "(the, pattern)" --src_tree_list trees.file --root | wc -l`


The render option will save each match as an image. If there are multiple patterns, numbers will be used to designate each pattern starting from 0.
If there are multiple matches, and underscore is used with a number for each match starting with 0. If I had two

` ete_search --pattern_tree_list "MyPatterns.txt" --tree_format 8 --src_tree_list "MyTargetTrees.txt" --render treematches.png `
