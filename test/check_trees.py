#!/usr/bin/env python

import re
import sys
from dendropy import Tree, TaxonNamespace, TreeList
from dendropy.calculate import treecompare

def main():
    treefile1 = sys.argv[1]
    treefile2 = sys.argv[2]

    treelist = TreeList()
    treelist.read(file=open(treefile1, 'rU'), schema="nexus")
    treelist.read(file=open(treefile2, 'rU'), schema="nexus")
        
    if treecompare.symmetric_difference(treelist.__getitem__(0), treelist.__getitem__(1)) == 0:
        print "trees are identical"
    else:
        print "trees are NOT identical"

if __name__ == '__main__':
    main()
