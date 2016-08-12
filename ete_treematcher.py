import sys
from logging import log
from argparse import ArgumentParser
from .common import src_tree_iterator
from ..phylo import PhyloTree
from treematcher import TreePattern, TreePatternCache

DESC=''

#ete3 treematcher --pattern "(hello, kk);" --pattern-format 8 --tree-format 8 --trees "(hello,(1,2,3)kk);" --quoted-node-names


def populate_args(extract_args_p):

    extract_args = extract_args_p.add_argument_group('TREEMATCHER OPTIONS')

    extract_args.add_argument("-p", dest="pattern",
                              nargs=1,
                              help="Supply a pattern you are going to search for.")
    extract_args.add_argument("--quoted-node-names", dest="quoted_node_names",
                              action="store_true",
                              help="True if using quotes to designate node names in the pattern. Otherwise False.")
    extract_args.add_argument("--tree-format", dest="tree_format",
                              type=int,
                              default=0,
                              help="A number 0-8 designating Newick format.")
    extract_args.add_argument("--maxhits", dest="maxhits",
                              nargs=1,
                              type=int,
                              default=0,
                              help="The number of matches to return. Default is 0 which returns all matches.")
    extract_args.add_argument("--cache", dest="cache",
                              action="store_true",
                              help="True if a cache is to be used. Otherwise False.")

def run(args):
    if args.src_trees is None:
        log.error('Please specify a tree to search (i.e. -t) ')
        sys.exit(-1)
    if args.pattern is None:
        log.error('Please specify a pattern to search for. (i.e. -p)')
        sys.exit(-1)
    pattern = TreePattern(args.pattern[0])

    if args.maxhits == 0:
        args.maxhits = None

    if args.verbosity > 2:
        print("Pattern is: ")
        print(pattern)

    if args.output:
        outputfile = open(args.output, 'w')

    for nw in src_tree_iterator(args):
        t = PhyloTree(nw, format=args.tree_format,
                     quoted_node_names=args.quoted_node_names)
        if args.cache:
            cache = TreePatternCache(t)
        else:
            cache = None
        matches = list(pattern.find_match(t, maxhits=args.maxhits, cache=cache))
        if args.output:
            outputfile.write(str(matches) + '\n')
        else:
            print(matches)

    if args.output:
        outputfile.close()



if __name__ == "__main__":
    parser = ArgumentParser()
    populate_args(parser)
    args = parser.parse_args(sys.argv[1:])
    run(args)



    '''
    #add in ete.py

     from . import (ete_split, ete_expand, ete_annotate, ete_ncbiquery, ete_view,
               ete_generate, ete_mod, ete_extract, ete_compare, ete_evol,
               ete_maptrees, ete_treematcher)


    # -treematcher-
    treematcher_args_p = subparser.add_parser("treematcher", parents=[source_args_p, main_args_p],
                                           description=ete_treematcher.DESC)
    treematcher_args_p.set_defaults(func=ete_treematcher.run)
    ete_treematcher.populate_args(treematcher_args_p)


    '''