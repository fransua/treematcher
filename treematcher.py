from itertools import permutations
from string import strip
import re
import sys

from ete3 import PhyloTree, Tree

class TreePattern(Tree):
    def __init__(self, *args, **kargs):
        kargs["format"] = 1
        Tree.__init__(self, *args, **kargs)
        for n in self.traverse():
            if n.name != "NoName":
                n.constraint = n.name.replace("{", "(").replace("}", ")").replace("@", "__target").replace("|", ",")
            else:
                n.constraint = None

    def constrain_match(self, __target, local_vars = None):
        if not self.constraint:
            return True
        
        if not local_vars:
            local_vars = {}
        local_vars.update({"__target": __target, "self": __target})
        try:
            st = eval(self.constraint, local_vars) if self.constraint else True
            #print __target
            st = bool(st)
        except ValueError: 
                raise ValueError("The following constraint expression did not return boolean result: %s BUT %s" %
                                 (self.constraint, st))

        return st
    
    def is_match(self, node, local_vars=None):
        # Check expected features
        status = self.constrain_match(node, local_vars)
        if status and self.children:
            #print "has children"
            if len(node.children) >= len(self.children):
                # Check all possible comparison between pattern children and
                # and tree node children.
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
        return self.get_ascii(show_internal=True, attributes=["constraint"])

    def find_match(self, tree, local_vars):
        for node in tree.traverse("preorder"):
            if self.is_match(node, local_vars=local_vars):
                return True, node
        return False, None


        
    
def length(txt):
    return len(txt)


def test():

    # custom_functions = {"length": length, "optional": optional}
    custom_functions = {"length": length}


    pattern = """
    (len{@.children} > 2 and @.name in {"hello"|"bye"}){length{@.name} < 3 or @.name == "pasa"} and @.dist >= 0.5;
    """

    pattern = TreePattern(pattern, format=8)

    print pattern
    '''
    tree = Tree("((hello,(1,2,3)kk)pasa:1, NODE);", format=1)
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree?:", pattern.find_match(tree, custom_functions)

    tree = Tree("((hello,(1,2,3)kk)pasa:0.4, NODE);", format=1)
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree?:", pattern.find_match(tree, custom_functions)

    tree = Tree("(hello,(1,2,3)kk)pasa:1;", format=1)
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree?:", pattern.find_match(tree, custom_functions)

    tree = Tree("((bye,(1,2,3)kk)none:1, NODE);", format=1)
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree?:", pattern.find_match(tree, custom_functions)
    '''
    print "looking at tree 1."
    tree = Tree("((kk,(1,2,3)bye)y:1, NODE);", format=1)
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree?:", pattern.find_match(tree, custom_functions)

    print "looking at tree 2"
    tree = Tree("(((1,2,3)bye)y:1, NODE);", format=1)
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree?:", pattern.find_match(tree, custom_functions)

    print "looking at tree 3."
    tree = Tree("(((1,2,3)bye,kk)y:1, NODE);", format=1)
    print tree.get_ascii(attributes=["name", "dist"])
    print "Pattern matches tree?:", pattern.find_match(tree, custom_functions)


if __name__ == "__main__":
    test()