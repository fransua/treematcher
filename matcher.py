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

def _parse_pattern(node, NHX_string):
    """ Reads node's extra data form its NHX string. NHX uses this
    format:  [&&NHX:prop1=value1:prop2=value2] """
    NHX_string = NHX_string.replace("[&&NHX:", "")
    NHX_string = NHX_string.replace("]", "")
    for field in NHX_string.split(","):
	ok = False
	for key, (op, cast) in OPERATOR.iteritems(): 
	    try:
		attr, value = map(strip, field.split(key))
	    except ValueError, e:
                pass
	    else:
		node.conditions.append([attr, op, value])
		ok = True
		break
	if not ok: 
	    raise ValueError("Unknown operator [%s]" %field)

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
	self.conditions = []	    
	if newick is not None:
	    nwParser._parse_extra_features = _parse_pattern
	    nwParser.read_newick(newick, root_node = self, format=format)

OPERATOR = {
    ">=": [OP.ge, float], 
    "<=": [OP.le, float],
    ">": [OP.gt, float],
    "<": [OP.lt, float],
    "==": [OP.eq, float],
    "~": [like, str],
    }

def _match(p, t):
    for attr, op, expected in p.conditions: 
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
    
    # Enter into topology comparison
    matches = []
    if status and len(P.children):
        if len(T.children) == len(P.children):
            # Check all possible comparison between patter children and
            # and tree node children.
            for candidate in permutations(P.children):
                sub_status = True
                for i, ch in enumerate(candidate): 
                    st = match(ch, T.children[i]) 
                    sub_status &= st
                status = sub_status
        else:
            status = False

    return status


pattern = "(([&&NHX:species==Mmu],([&&NHX:species==Hsa],[&&NHX:species==Hsa])),[&&NHX:species~Pt.]);"
nw = "((Mms01, (Hsa01, Ptr01)), ((Mmu01, (Hsa01, Hsa02)),Ptr02), ((Mmu01, (Hsa01, Hsa02)),Pts02) );"


P = TreePattern(pattern)

print P.write(features=[])

T = PhyloTree(nw)

for n in T.traverse():
    if match(P, n):
	print "MATCHING NODE"
	print n

CountA = CountLeaves(lambda x: x.name == "A")
isLeaf = CheckLeaf()
pattern = "(any[&&NHX:isLeaf{}==1],moreA[&&NHX:CountA{}>=2]);"

P = TreePattern(pattern)
print P
for n in P.traverse():
    print n.conditions
print len(P.children)
T = PhyloTree("(((A,B),(C,D)),U);")
#print T
T.populate(200, names_library="ABC", reuse_names=True)
print T
for n in T.traverse():
    if match(P, n):
	print "MATCHING NODE"
	print n
        raw_input()


