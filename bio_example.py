

"""
 The trees used in the following examples were taken from
 Morgan, Claire C., Christopher J. Creevey, and Mary J. O'Connell.
 Mitochondrial data are not suitable for resolving placental mammal phylogeny.
 Mammalian Genome 25.11-12 (2014): 636-647.
"""
from ast import literal_eval
from ete3 import NCBITaxa, PhyloTree
from treematcher import TreePattern

ncbi = NCBITaxa()

with open("tree_examples.txt", 'r') as fh:
    tree_dict = literal_eval(fh.read())  # read in example trees

#########################################################################
# Example 1
# Suppose we want to double check that all of the trees listed under "Insectivora"
# have leaves which are classified in Insectivora according to the NCBI Taxonomy database.
# We expect one outgroup, Balaena_mysticetus, and we want to check if it is
# placed correctly in the tree.
border = "\n" + "#" * 110 + "\n"
print(border)
print("Example 1: Check that leaves are classified in Insectivora according to the NCBI Taxonomy database.")
print(" We expect one outgroup, Balaena_mysticetus, and we want to check that it diverges from root. ")
print(border)
#########################################################################

Insectivora_trees = tree_dict["Insectivora"]  # Access Insectivora dictionary of trees
Insectivora_trees = [PhyloTree(tree) for tree in Insectivora_trees]  # creates tree pattern instances

# Find all leaves in tree that are not Insectivora
pattern1 = TreePattern(""" '  @.is_leaf() and "Insectivora" not in @.named_lineage ';""")

# Check that outgroup is first to diverge from root
pattern2 = TreePattern("""  ( 'Balaena_mysticetus' ) ' @.is_root() ';""")

count=1

for tree in Insectivora_trees:
    # set species
    tree.set_species_naming_function(
        lambda n: (ncbi.get_name_translator([n.name.replace("_", " ")])).get(n.name.replace("_", " "))[
            0] if "_" in n.name else '')

    # access taxonomy information
    tree.annotate_ncbi_taxa()

    # Find non-Insectivora, check if outgroup, look at placement
    for match in list(pattern1.find_match(tree, maxhits=None)):
        if len(match)==1:
            print("All leaves in tree {} are Insectivora".format(count))
        if match.name=="Balaena_mysticetus":  # outgroup
            print("outgroup diverges from root? {}  ".format(len(list(pattern2.find_match(tree))) == 1))
            # True if outgroup is child of root
        else:
            print("Not Insectivora: {} , tree number: {} ".format(match.name, count))
    count+=1

# All trees return the correct outgroup
# but only the first tree has the outgroup diverging directly from the root

#########################################################################
# Example 2
# Suppose a different arrangement of orders was published in a new study, such as
# -- /-Insectivora/-Chiroptera/-Cetartiodactyla/-Perissodactyla/-

# and we want to check if this arrangement is present within the trees
# listed under Laurasiatheria
print(border)
print("Example 2: Looking for the presence of nodes in the arrangement")
print("-- /-Insectivora/-Chiroptera/-Cetartiodactyla/-Perissodactyla/-Carnivora")
print(border)
#########################################################################

Laurasaitheria_trees = tree_dict["Laurasaitheria"]
Laurasaitheria_trees = [PhyloTree(tree) for tree in Laurasaitheria_trees]

pattern1 = TreePattern(""" '  @.sci_name == "Insectivora" ';""")
pattern2 = TreePattern(""" '  @.sci_name == "Chiroptera" ';""")
pattern3 = TreePattern(""" '  @.sci_name == "Cetartiodactyla" ';""")
pattern4 = TreePattern(""" '  @.sci_name == "Perissodactyla" ';""")
pattern5 = TreePattern(""" '  @.sci_name == "Carnivora" ';""")
count = 1
for tree in Laurasaitheria_trees:
    # set species attribute
    tree.set_species_naming_function(
        lambda n: (ncbi.get_name_translator([n.name.replace("_", " ")])).get(n.name.replace("_", " "))[
            0] if "_" in n.name else '')

    # access taxonomy information from NCBI database
    tree.annotate_ncbi_taxa()

    match1 = list(pattern1.find_match(tree, maxhits=1, target_traversal="levelorder"))
    match2 = list(pattern2.find_match(tree, maxhits=1, target_traversal="levelorder"))
    match3 = list(pattern3.find_match(tree, maxhits=1,target_traversal="levelorder"))
    match4 = list(pattern4.find_match(tree, maxhits=1, target_traversal="levelorder"))
    match5 = list(pattern5.find_match(tree, maxhits=1, target_traversal="levelorder"))

    results = len(match1) + len(match2) + len(match3) + len(match4) + len(match5)

    if results>=2:
        try:
            if match1:  # if Insectivora exists in tree, make it the new subtree to search
                tree = match1[0].up

            if match2:  # if Cetartiodactyla exists in tree
                # check that it exists in the current subtree
                subtree = list(pattern2.find_match(tree, maxhits=1, target_traversal="levelorder"))[0]
                if subtree != tree:
                    tree = subtree.up
                else:
                    raise IndexError

            if match3:  # if Perissodactyla  exists in tree
                # check that it exists in the current subtree
                subtree = list(pattern3.find_match(tree, maxhits=1, target_traversal="levelorder"))[0]
                if subtree!=tree:
                    tree = subtree.up
                else:
                    raise IndexError

            if match4: # if carnivora exists, check that it exists in the subtree
                subtree = list(pattern4.find_match(tree, maxhits=1, target_traversal="levelorder"))[0]
                if subtree != tree:
                    tree = subtree.up
                else:
                    raise IndexError

            if match5:  # if carnivora exists, check that it exists in the subtree
                subtree = list(pattern4.find_match(tree, maxhits=1, target_traversal="levelorder"))[0]
                if subtree == tree:
                    raise IndexError

            print(" {} orders that are present occur in the correct arrangement for tree {} ".format(results, count))

        except IndexError:
            print(" {} orders not the correct arrangement for tree {}".format(results, count))  # no solution for order on subtree
    else:
        print("Only one Order found in tree {}".format(count))
    count+=1
# None of the trees contain all of the groups.
# Trees 2,3,10 and 12 have the correct for the orders that are present