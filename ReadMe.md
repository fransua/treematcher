# TreeMatcher: A new tool for creating regular-expression-like queries on trees

In mathematics, a standard way of representing graphical trees with edge lengths is the Newick format which uses a specific syntax (such as parentheses and commas). The TreeMatcher module defines a tree pattern by extending the Newick format to include rules and filters with a regular-expression-like vocabulary. These patterns are then searched for using a tree traversal algorithm.

### Syntax

The treematcher syntax consists of all of the attributes and functions available to a TreeNode written inside a Newick structure. The pattern is usually begins as a string surrounded by three double quotes at each end which is used to create TreePattern instance. To differentiate the parentheses of a set, tuple, or function in your pattern verses the parentheses of Newick format, use quoted node names.


Examples


|  syntax       						| meaning       						|  Comments																|
| ----------------------------------	|:----------------------------------:	|:---------------------------------------------------------------------:|
| @	            						| tree node								| Use @.attribute to access an attribute, function(@) to access function|
| @.name == "sample1"					| equivalent to @.name == "sample1" 	| Looking for multiple names, use list: @.name in ("sample1","sample2") |
| @.dist >= 0.5     					| branch length no less than 0.5		| 	Use any of the following: <, <=, ==, >=								|
| @.support > 0.9	            		| Has a support value greater than 0.90	| 																		|
| @.species=="Homo sapiens"	    		| Homo sapiens is species of node		| See set_species_naming_function()	for details							|
| 9606 in @.lineage	            		| Homo sapiens in @.lineage				| Find NCBI taxonomy ID or the full scientific name	in a node's lineage	|
| @.sci_name == "Euarchontoglires" 		| scientific name is Euarchontoglires	| See annotate_ncbi_taxa() function for details 						|
| @.rank=="subfamily" 					| node is ranked at the subfamily level	| See annotate_ncbi_taxa() function for details							|
| @.taxid == 207598						|20758 is the taxid	of the node			| See annotate_ncbi_taxa() function for details							|
| "Pan_troglodytes_1" in @.leaves		| Pan_troglodytes_1 is descendant leaf	| Find the leaf name within a list of leaf names						|
| len(@.children)						| number of children					| use quotes to differentiate parentheses from Newick Structure			|
| "Homo sapiens" in @.contains_species	| find species contained in a leaf		| custom attribute														|
| @.size > 5							| number of descendants/size of tree	| custom attribute 														|
| "H_saps_1" in @.leaves				| find name H_saps_1 in leaves			| custom attribute														|


### How to use Newick format to access the structure of a tree
Use Newick format to access the tree structure. Call the TreePattern function to create an instance of a tree pattern.


Example 1: Find a tree where sample_1 and sample_2 are siblings.
```
pattern1 = """ sample_1, sample_2 ; """
pattern1 = TreePattern(pattern1, format=8, quoted_node_names=True)
```

Example 2: Find a tree where sample_1 and sample_2 are children of the parent sample_0.
```
pattern2 = """ (sample_1, sample_2) sample_0 ; """
pattern2 = TreePattern(pattern0, format=8, quoted_node_names=True)
```

### To Run
To run, use the find_match function.
Example1: Find the node that matches pattern1 in the following tree:

tree = Tree("(sample_1,(sample_1,sample_2)sample_0)sample_0:1;", format=1)

```
pattern1.find_match(tree, None)
```

Example2: Find the total number of pattern2 matches in the same tree
``` len(pattern2.find_match(tree, None, maxhits=None) ```


### custom functions
Write your own functions and provide them as local variables to the treematcher program. You can alter the function names using a dictionary.
```
pattern3 = """( 'contains(@, ("Chimp_2", "Chimp_3"))' , 'num_species(@, 2) and num_leaves(@,2)' ); """

pattern3 = TreePattern(pattern3, format=8, quoted_node_names=True,
                      functions={'contains': contains,
                                 "num_species": number_of_species,
                                 "num_leaves": number_of_leaves})

pattern3.find_match(tree, None, maxhits=None)
```
