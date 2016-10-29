#!/usr/bin/env python

import re
import sys
import xmltodict

def main():
    treefile1 = sys.argv[1]
    treefile2 = sys.argv[2]

    f = open(treefile1, 'r')
    xml = f.read()
    f.close()
    xmldict1 = xmltodict.parse(xml, force_list=('tree',))

    f = open(treefile2, 'r')
    xml = f.read()
    f.close()
    xmldict2 = xmltodict.parse(xml, force_list=('tree',))
    
    tree1 = get_tree_from_nexml(xmldict1)
    tree2 = get_tree_from_nexml(xmldict2)
    same_tree(tree1, tree2)

def same_tree(tree1, tree2):
    print tree1
    print tree2
    if len(tree1['children']) == 0 and len(tree2['children']) == 0:
        # if both trees are otus, they're equal
        return True
    elif len(tree1['children']) == 0:
        # if one is an otu and the other isn't, they're not.
        return False
    elif len(tree2['children']) == 0:
        return False
    else:
        #augh this might be hard: for each child, see if one of the other tree's children matches.
        #if 
        print "hi"
# 1. If both trees are empty then return 1.
# 2. Else If both trees are non -empty
#      (a) Check data of the root nodes (tree1->data ==  tree2->data)
#      (b) Check left subtrees recursively  i.e., call sameTree( 
#           tree1->left_subtree, tree2->left_subtree)
#      (c) Check right subtrees recursively  i.e., call sameTree( 
#           tree1->right_subtree, tree2->right_subtree)
#      (d) If a,b and c are true then return 1.
# 3  Else return 0 (one is empty and other is not)

    
def get_tree_from_nexml(nexmldict):
    trees = nexmldict['nex:nexml']['trees']
    currtree = trees['tree'][0]
    otus = {}
    nodes = []
    edges = []
    for node in currtree['node']:
        if '@otu' in node:
            otus[node['@otu']] = node['@id']
            nodes.append(node['@id'])
        else:
            nodes.append(node['@id'])
    for edge in currtree['edge']:
        edges.append({'source':edge['@source'],'target':edge['@target']})

    nodedict = {}
    for k in otus.keys():
        nodedict[otus[k]] = []
    # find the root: it's the node that is not the target of any edge
    root = None
    for node in nodes:
        targetnum = 0
        sourcenum = 0
        if node not in nodedict:
            nodedict[node] = []
        for edge in edges:
            if node == edge['source']:
                targetnum += 1
                edgenex = str(edge['target'])
                nodedict[node].append(edgenex)
            if node == edge['target']:
                sourcenum += 1
        if sourcenum == 0:
            root = {'node':node,'children':nodedict[node]}
    final_tree = replace_children(nodedict, root)
    return final_tree

def replace_children(nodedict, node):
    new_children = []
    if len(node['children']) == 0:
        return node
    else:
        for child in node['children']:
            new_node = {'node':child, 'children':nodedict[child]}
            new_children.append(new_node)
            replace_children(nodedict, new_node)
    node['children'] = new_children
    return node

# def write_tree_as_string(tree):
    # wa

if __name__ == '__main__':
    main()
