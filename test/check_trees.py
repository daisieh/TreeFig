#!/usr/bin/env python

import re
import sys

def main():
    treefile1 = sys.argv[1]
    treefile2 = sys.argv[2]

    f = open(treefile1, 'r')
    nexus1 = f.read()
    f.close()

    f = open(treefile2, 'r')
    nexus2 = f.read()
    f.close()
    
    tree1 = get_newick(nexus1)
    tree2 = get_newick(nexus2)
    print same_tree(tree1, tree2)

def get_newick(str):
    treematcher = re.search('Tree\s+.+?=\s+(.+?);', re.sub('\n',' ', str), re.MULTILINE)
    if treematcher is not None:
        return treematcher.group(1)
    return ''
    
# sameTree(tree1, tree2)
# 1. If both trees are empty then return 1.
# 2. Else If both trees are non -empty
#      (a) Check data of the root nodes (tree1->data ==  tree2->data)
#      (b) Check left subtrees recursively  i.e., call sameTree( 
#           tree1->left_subtree, tree2->left_subtree)
#      (c) Check right subtrees recursively  i.e., call sameTree( 
#           tree1->right_subtree, tree2->right_subtree)
#      (d) If a,b and c are true then return 1.
# 3  Else return 0 (one is empty and other is not)
def same_tree(tree1, tree2):
    if tree1 == tree2:
        return True
    treedata1 = []
    treedata2 = []
    pattern = re.compile('\((.*)\)')
    treedatamatcher = pattern.match(tree1)
    if treedatamatcher is not None:
        treedata1.append(treedatamatcher.group(1))
    else:
        nodes = re.split(',',tree1)
        treedata1.extend(nodes)
    treedatamatcher = pattern.match(tree2)
    if treedatamatcher is not None:
        treedata2.append(treedatamatcher.group(1))
    else:
        nodes = re.split(',',tree2)
        treedata2.extend(nodes)

    if treedata1 == treedata2:
        print "%s equals %s" % (treedata1, treedata2)
        return True
    if len(treedata1) != len(treedata2):
        return False
    
    
    return same_tree(treedata1[0], treedata2[0])    
    return False

if __name__ == '__main__':
    main()
