from itertools import permutations
from string import strip
import re
from ete2 import PhyloTree, Tree
from ete2.parser import newick as nwParser
import operator as OP

def like(pattern,string):
    return re.search(string, pattern)

class CountLeaves(object):
    def __init__(self, filter_fn=None):
        self.filter_fn = filter_fn

    def run(self, node):
        i = 0
        for n in node.iter_leaves():
            if self.filter_fn(n):
                i+=1
        return i

class CheckLeaf(object):
    def run(self, node):
        return node.is_leaf()

isLeaf = CheckLeaf()

def _parse_pattern(node, NHX_string):
    """ Reads node's extra data form its NHX string. NHX uses this
    format:  [&&NHX:prop1=value1:prop2=value2] """
    NHX_string = NHX_string.replace("[&&NHX:", "")
    NHX_string = NHX_string.replace("]", "")
    is_leaf_check = False
    for field in NHX_string.split(":"):
        if field == "!":
            is_leaf_check = True
            continue
           
	ok = False
	for key, op in OPERATOR.iteritems(): 
	    try:
		attr, value = map(strip, field.split(key))
	    except ValueError, e:
                pass
	    else:
		node.constrain.append([attr, key, value])
		ok = True
		break
	if not ok: 
	    raise ValueError("Unknown operator [%s]" %field)

    if is_leaf_check:
        node.constrain.append(["isLeaf{}", "==", True])


class TreePattern(Tree):
    def __init__(self, newick=None, format=0):
	self._children = []
	self._up = None
	self._dist = 1.0
	self._support = 1.0
	self.name = ""
	self.features = set([])
	self.collapsed = False
	# New for pattern type
	self.constrain = []	    
	if newick is not None:
	    nwParser._parse_extra_features = _parse_pattern
	    nwParser.read_newick(newick, root_node = self, format=format)


    def __str__(self):
        return self.get_ascii(show_internal="True")

    def _asciiArt(self, char1='-', show_internal=True, compact=False):
        """
        Returns the ASCII representation of the tree. Code taken from the
        PyCogent GPL project.
        """
        constrain = ", ".join(map(lambda x: "%s %s %s" %(x[0],x[1],x[2]), 
                                 [(attr, op, value) for attr, op, value in self.constrain]))
        name_txt = self.name + " (%s)" %constrain

        LEN = 5
        PAD = ' ' * LEN
        PA = ' ' * (LEN-1)
        if not self.is_leaf():
            mids = []
            result = []
            for c in self.children:
                if c is self.children[0]:
                    char2 = '/'
                elif c is self.children[-1]:
                    char2 = '\\'
                else:
                    char2 = '-'
                (clines, mid) = c._asciiArt(char2, show_internal, compact)
                mids.append(mid+len(result))
                result.extend(clines)
                if not compact:
                    result.append('')
            if not compact:
                result.pop()
            (lo, hi, end) = (mids[0], mids[-1], len(result))
            prefixes = [PAD] * (lo+1) + [PA+'|'] * (hi-lo-1) + [PAD] * (end-hi)
            mid = (lo + hi) / 2
            prefixes[mid] = char1 + '-'*(LEN-2) + prefixes[mid][-1]
            result = [p+l for (p,l) in zip(prefixes, result)]
            if show_internal:
                stem = result[mid]
                result[mid] = stem[0] + name_txt + stem[len(name_txt)+1:]
            return (result, mid)
        else:
            return ([char1 + '-' + name_txt], 0)

OPERATOR = {
    ">=": OP.ge,
    "<=": OP.le,
    ">": OP.gt,
    "<": OP.lt,
    "==": OP.eq,
    "~": like,
    }

def _match(p, t):
    for attr, keyop, expected in p.constrain: 
        op = OPERATOR[keyop]
        if attr.endswith("{}"):
            fn = globals()[attr.strip("{}")]
            value = fn.run(t)
        else:
            value = getattr(t,attr)
        expected = type(value)(expected)
        if not op(value, expected):
            return False
    return True

def match(P, T):
    # Check expected features
    status = _match(P, T)
    #print P, T
    #print "constrains", status
    # Enter into topology comparison
    matches = []
    if status and len(P.children):
        #print "has children"
        if len(T.children) == len(P.children):
            # Check all possible comparison between patter children and
            # and tree node children.
            for candidate in permutations(P.children):
                sub_status = True
                for i, ch in enumerate(candidate): 
                    st = match(ch, T.children[i]) 
                    sub_status &= st
                status = sub_status
                if status:
                    break
        else:
            #print "!= length"
            status = False

    return status


CountA = CountLeaves(lambda x: x.name == "A")


pattern = "(([&&NHX:!:species==Mmu],([&&NHX:!:species==Hsa],[&&NHX:!:species==Hsa])),[&&NHX:!:species~Pt.]);"
nw = "((Mms01, (Hsa01, Ptr01)), ((Mmu01, (Hsa01, Hsa02)),Ptr02), ((Mmu01, (Hsa01, Hsa02)),Pts02) );"


P = TreePattern(pattern)

print P.write(features=[])

T = PhyloTree(nw)
print T, P
for n in T.traverse():
    if match(P, n):
	print "MATCHING NODE"
	print n
raw_input("continue")

pattern = "(any, moreA[&&NHX:CountA{}>=1]);"

P = TreePattern(pattern)
print P
T = PhyloTree("(((A,B),(C,D)),U);")
#print T
#T.populate(200, names_library="ABC", reuse_names=True)
print T
for n in T.traverse():
    if match(P, n):
	print "MATCHING NODE"
	print n
        raw_input()


