# TreeMatcher: A new tool for creating regular-expression-like queries on trees

In mathematics, a standard way of representing graphical trees with edge lengths is the Newick format which uses a specific syntax (such as parentheses and commas). The TreeMatcher module defines a tree pattern by extending the Newick format to include rules and filters with a regular-expression-like vocabulary. These patterns are then searched for using a tree traversal algorithm.

### Syntax

ETE uses its own formalism to represent phylogenetic trees. A pattern can be written by accessing the attributes and functions available to an ETE tree, using Python code, or through custom functions and the pattern structure is created through the Newick format.
The pattern is usually begins as a string surrounded by three double quotes at each end which is then used to create TreePattern instance. To differentiate the parentheses of a set, tuple, or function in your pattern verses the parentheses of Newick format, use quoted node names.


Examples


|  syntax       						| meaning       						|  Comments																|
| ----------------------------------	|:----------------------------------:	|:---------------------------------------------------------------------:|
| @	            						| node, default for nodes left blank	| Use @.attribute to access an attribute, function(@) to access function|
| "sample1"								| equivalent to @.name == "sample1" 	| Looking for multiple names, use list: @.name in ("sample1","sample2") |
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
Newick format access a tree structure. Use TreePattern (or its superclass tree) attributes and functions for the constraints on each node.
Call the TreePattern function to create an instance of a tree pattern.


Example 1: Find a tree where sample_1 and sample_2 are siblings.
```
pattern1 = """ sample_1, sample_2 ; """
pattern1 = TreePattern(pattern1)
```

Example 2: Find a tree where sample_1 and sample_2 are children of the parent sample_0.  Note that the format type is set to 1 as the default. You can access other Newick format types by setting the format argument.
```
pattern2 = """ (sample_1, sample_2) sample_0 ; """
pattern2 = TreePattern(pattern0, format=8)
```

### To Run
To run, use the find_match function.


Example1: For the following tree, find the node that matches pattern1.

tree = Tree("(sample_1,(sample_1,sample_2)sample_0)sample_0:1;", format = 8)

```
pattern1.find_match(tree, None)
```

Example2: Find the total number of pattern2 matches in the same tree
``` len(pattern2.find_match(tree, None, maxhits=None) ```


### custom functions
Write your own functions and provide them as local variables to the treematcher program. You can alter the function names using a dictionary. For example, suppose you have written two functions: num_of_species and num_of_leaves. Access these functions as follows.
```
pattern3 = """ 'sn(@, 2) and ln(@,2)' ; """

pattern3 = TreePattern(pattern3, format=8, quoted_node_names=True,
                      functions={"sn": number_of_species,
                                 "ln": number_of_leaves})

pattern3.find_match(tree, None, maxhits=None)
```
