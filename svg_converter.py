#!/usr/bin/env python

import re
import xmltodict
import sys
import json
import os
import math
from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier, parse_path

from dendropy import Tree, TaxonNamespace

total_width = 0
total_height = 0
scale_width = 1.0
scale_height = 1.0
max_x = 0
max_y = 0
min_x = 0
min_y = 0
otu_level = 0
root_level = 0

points = []

def main():
    global points
    filename = sys.argv[1]
    # convert input_image -threshold 50% raw.pbm
    # potrace -s -k 0.8 -W 10 -H 10 -o output.svg raw.pbm
    file_name, extension = os.path.splitext(sys.argv[1])
    outfile = file_name
    if extension != '.svg':
        print "can't open this file"
        return
    f = open(filename, 'r')
    xmldict = ''
    xml = f.read()
    f.close()
    xmldict = xmltodict.parse(xml, force_list=('path',))
    xmldict = xmldict['svg']
    global total_width, total_height, scale_width, scale_height
    s = re.sub('[a-zA-Z]+', '', xmldict['@width'])
    total_width = float(s)
    s = re.sub('[a-zA-Z]+', '', xmldict['@height'])
    total_height = float(s)
    xmlpaths = []
    style = {}
    transform = ''
    paths = [] 
    circles = []

    # parse original paths
    if 'path' in xmldict:
        xmlpaths = xmldict['path']
    elif 'g' in xmldict:
        xmlpaths = xmldict['g']['path']
        del xmldict['g']['path']
        if '@transform' in xmldict['g']:
            transform = xmldict['g']['@transform']
            if transform != '':
                parse_transform(transform)
    raw_polygons = []
    for path in xmlpaths:
        polygon = path_to_polygon(path['@d'])
        raw_polygons.append(polygon)

    global max_x, max_y, min_x, min_y
    global scale_width, scale_height
    max_x = abs(max_x * scale_width)
    max_y = abs(max_y * scale_height)
    treepaths = []
    otherpaths = []
    boxpaths = []
    for polygon in raw_polygons:
        polygon = scale(polygon, scale_width, scale_height)
        polygon = translate(polygon, 0, max_y)
        box = bounding_box(polygon)
        if threshold_area(box, 0.6):
            path = {}
            path['@d'] = nodes_to_path(polygon)
            path['@style'] = "fill:#CCCCCC; stroke:#999999; stroke-width:1"
            otherpaths.append(path) 
            treepaths.append(polygon)
        else:
            path = {}
            path['@d'] = nodes_to_path(polygon)
            path['@style'] = "fill:#EEEEEE; stroke:#EEEEEE; stroke-width:1"
            paths.append(path) 
            otherpaths.append(path)
            boxpaths.append(box)

    boxpaths = process_boxpaths(boxpaths)
    i = 0
    for box in boxpaths:
        path = {}
        path['@d'] = nodes_to_path(box)
        val = i%256
        i = i+10
        path['@style'] = "fill:#00DD%x; stroke:#EEEEEE; stroke-width:1" % val
        otherpaths.append(path)

    segments = []
    rawtreepaths = []
    rawpolygons = []
    otus = []
    for treepath in treepaths:
        node_distance = find_node_distance(treepath)
        scaling = 50.0/node_distance
        print "node distance is %d, scaling is %d" % (node_distance, scaling)
        for node in treepath:
            node[1] = node[1] * scaling
        straighten_horizontal_lines(treepath)
        straighten_polygon(treepath)
        # pointify_tips(treepath)
        otus.extend(find_otus(treepath))

        # for node in treepath:
        #     node[1] = node[1] / scaling

        # this path is for the cleaned-up lines
        path = {}
        path['@d'] = nodes_to_path(treepath)
        path['@name'] = "cleaned path"
        path['@style'] = "fill:none; stroke:#FF0000; stroke-width:1"
        rawtreepaths.append(path) 

        segments.extend(polygon_to_lines(treepath))
    # generate raw svg first-pass, in case something fails during tree building:
    print "otus " + str(otus)
    circles.extend(nodes_to_circles(otus))
    lines = []
    lines.extend(segments_to_lines(segments, 'blue', 4))
    svgdict = {}
    svgdict['svg'] = {}
    svgdict['svg']['@width'] = xmldict['@width']
    svgdict['svg']['@height'] = xmldict['@height']
    svgdict['svg']['g'] = []
#     svgdict['svg']['g'].extend([{'line':lines}])
    svgdict['svg']['g'].extend([{'circle':circles}])
    svgdict['svg']['g'].extend([{'path':otherpaths}])
    svgdict['svg']['g'].extend([{'path':rawtreepaths}])


    outf = open(outfile+'_raw.svg','w')
    outf.write(xmltodict.unparse(svgdict, pretty=True))
    outf.close()
    
    print "printed raw svg"
    
    # make the raw tree:
    (nodes, edges, otus) = make_tree(segments)
    newick_tree = tree_to_nexus(otus, nodes, edges) + ';'
    dendro_tree = Tree.get(data=newick_tree, schema="newick")
    
    # generate nexus:
    dendro_tree.write(path=(outfile+'.nex'), schema="nexus")
    
    # generate nexml:
    dendro_tree.write(path=(outfile+'.xml'), schema="nexml")
            
    # process possible text labels:
    (otu_label_dict, branch_label_dict, other_label_dict) = define_text_boxes(boxpaths, nodes, segments, otus)
    
    for key in otu_label_dict:
        box = otu_label_dict[key]
        path = {}
        path['@d'] = nodes_to_path(box)
        path['@style'] = "fill:none; stroke:#3333DD; stroke-width:1"
        otherpaths.append(path)

    for key in branch_label_dict:
        box = branch_label_dict[key]
        path = {}
        path['@d'] = nodes_to_path(box)
        path['@style'] = "fill:none; stroke:#33DD33; stroke-width:1"
        otherpaths.append(path)

    for key in other_label_dict:
        boxes = other_label_dict[key]
        for box in boxes:
            path = {}
            path['@d'] = nodes_to_path(box)
            path['@style'] = "fill:#FFFF00; stroke:#33DD33; stroke-width:1"
            otherpaths.append(path)

        
    
    otudict = {}
    index = 1
    for otu in otus:
        otudict[str(otu)] = 'otu%d' % index
        index = index+1

    textdict = []
    # for each otu, add a text label:
    for k in otu_label_dict.keys():
        x = otu_label_dict[k][0][0]
        y = otu_label_dict[k][2][1]
        textnode = {'@x':x, '@y':y,'@fill':'black','#text':otudict[k]}
        textdict.append(textnode)
            
    # generate svg:
    nodes.extend(otus)
    circles.extend(nodes_to_circles(nodes))
    lines.extend(segments_to_lines(edges, "green", 3))
    # print lines
    svgdict = {}
    svgdict['svg'] = {}
    svgdict['svg']['@width'] = xmldict['@width']
    svgdict['svg']['@height'] = xmldict['@height']
    svgdict['svg']['g'] = [{'path':otherpaths}]
    svgdict['svg']['g'].extend([{'path':rawtreepaths}])
    svgdict['svg']['g'].extend([{'line':lines, 'circle':circles},{'text':textdict}])
    

#     outf = open(outfile+'_raw.svg','w')
#     outf.write(xmltodict.unparse(svgdict, pretty=True))
#     outf.close()
#     print "reprinted raw svg"


def find_node_distance(treepath):
    otus = find_otus(treepath)
    otus.sort(key=lambda r: r[1])
    print str(otus)
    shortest_distance = otus[1][1] - otus[0][1]
    current_node = otus.pop(0)
    while len(otus) > 0:
        if otus[0][1] - current_node[1] < shortest_distance:
            shortest_distance = otus[0][1] - current_node[1]
        current_node = otus.pop(0)
    return shortest_distance


def tree_to_nexus(otus, nodes, edges):
    otudict = {}
    nodedict = {}
    for i in range(len(otus)):
        nodedict[str(otus[i])] = str('otu%d' % i)
    # find the root: it's the node that is not the y2 of any edge
    root = None
    for node in nodes:
        targetnum = 0
        sourcenum = 0
        if str(node) not in nodedict:
            nodedict[str(node)] = []
        for edge in edges:
            if edge[1] == node[1] and edge[0] == node[0]:
                targetnum += 1
                edgenex = str([edge[2],edge[3]])
                nodedict[str(node)].append(edgenex)
            if edge[3] == node[1] and edge[2] == node[0]:
                sourcenum += 1
        if sourcenum == 0:
            root = node
    return replace_nodes(nodedict, str(root))

def replace_nodes(nodedict, newick):
    # find the nodes in newick:
    nodematcher = re.findall('\[\d+, \d+\]',newick)
    if len(nodematcher) is 0:
        return newick
    else:
        for node in nodematcher:
            # replace the node with the children of the node
            child = str(nodedict[node])
            if '[' in child:
                child = '(' + str(', '.join(nodedict[node])) + ')'
            newick = newick.replace(node, child)                
    return replace_nodes(nodedict, newick)

def scale(polygon, scale_width, scale_height):
    for node in polygon:
        node[0] = int(float(node[0]) * scale_width)
        node[1] = int(float(node[1]) * scale_height)
    return polygon
    
def translate(polygon, x, y):
    for node in polygon:
        node[0] = int(float(node[0]) + x)
        node[1] = int(float(node[1]) + y)
    return polygon

def parse_transform(transform):
    global scale_width, scale_height
    scalematcher = re.search('scale\(([0-9\-\.]+),*(.*?)\)', transform)
    if scalematcher is not None:
        scale_width = float(scalematcher.group(1))
        if scalematcher.group(2) is not None:
            scale_height = float(scalematcher.group(2))


def process_boxpaths(boxpaths):
    resultpaths = []
    # sort boxes by y:
    boxpaths.sort(key=lambda r: r[0][1])
    working_boxes = []
    # print "starting %s,%s" % (str(boxpaths[0]),str(boxpaths[1]))
    while len(boxpaths) > 0:
        if len(boxpaths) == 1:
            resultpaths.append(boxpaths.pop(0))
            boxpaths = working_boxes
            working_boxes = []
            # print "  finished a box %s" % str(resultpaths)
            if len(boxpaths) == 0:
                return resultpaths
        current_box = boxpaths.pop(0)
        comparison_box = boxpaths.pop(0)
        # print "  comparing boxes %s to %s" % (str(current_box),str(comparison_box))
        if compare_line_segments([current_box[0][1],current_box[1][1]],[comparison_box[0][1],comparison_box[1][1]]):
            boxpaths.insert(0, merge_boxes(current_box, comparison_box))
            # print "    %d boxes left" % len(boxpaths)
        else:
            boxpaths.insert(0, current_box)
            # print "    %d boxes left" % len(boxpaths)
            working_boxes.append(comparison_box)
        # break

    # print "boxpaths %s" % str(boxpaths)
    return resultpaths


def merge_boxes(box1, box2):
    # print "  merging boxes %s to %s" % (str(box1), str(box2))
    concat_points = box1 + box2
    concat_points.sort(key=lambda r: r[1])
    min_y = concat_points[0][1]
    max_y = concat_points[len(concat_points)-1][1]
    concat_points.sort(key=lambda r: r[0])
    min_x = concat_points[0][0]
    max_x = concat_points[len(concat_points)-1][0]
    return [[min_x,min_y],[min_x,max_y],[max_x,max_y],[max_x,min_y]]



def compare_line_segments(range1, range2):
    if range1[0] <= range2[0] <= range1[1]:
        return True
    if range1[0] <= range2[1] <= range1[1]:
        return True
    return False

def make_tree(segments):
    vert_lines = []
    horiz_lines = []
    for seg in segments:
        if seg[0] == seg[2]:
            vert_lines.append(seg)
        elif seg[1] == seg[3]:
            horiz_lines.append(seg)
    
    # create a dictionary of nodes: each node has vertical endpoints and is indexed by its x value
    node_dict = {}
    sorted_verts = sort_lines(vert_lines, 0, 1)

    for bin in sorted_verts:
        if bin[0][0] not in node_dict:
            node_dict[bin[0][0]] = []
        for line in bin:
            x = line[0]
            y1 = line[1]
            y2 = line[3]
            node_dict[x].append([y1,y2])
    
    # okay, now we know what the nodes are. Match up the edges.
    edges = set()
    otus = set()
    nodes = set()
    for line in horiz_lines:
        x1 = line[0]
        x2 = line[2]
        y1 = line[1]
        y2 = line[3]
        # we want to make the y1 equal to the y1 of the node
        # look for the node that this x1 is in:
        if x1 in node_dict:
            for node in node_dict[x1]:
                if y1 >= node[0] and y1 <= node[1]:
                    y1 = node[0]
            if x2 not in node_dict:
                otus.add('%d %d' % (x2, y2))
            else:
                for node in node_dict[x2]:
                    if y2 >= node[0] and y2 <= node[1]:
                        y2 = node[0]
            nodes.add('%d %d' % (x1, y1))
            nodes.add('%d %d' % (x2, y2))
            edges.add('%d %d %d %d' % (x1, y1, x2, y2))

    final_nodes = []
    for node in nodes:
        coords = re.split(' ',node)
        final_nodes.append([int(coords[0]), int(coords[1])])
    final_edges = []
    for edge in edges:
        coords = re.split(' ',edge)
        final_edges.append([int(coords[0]), int(coords[1]), int(coords[2]), int(coords[3])])
    final_otus = []
    for otu in otus:
        coords = re.split(' ',otu)
        final_otus.append([int(coords[0]), int(coords[1])])
    final_otus.sort(cmp=lambda x,y: cmp(x[1], y[1]))
    return (final_nodes, final_edges, final_otus)

def segments_to_lines(segments, color, width):
    lines = []
    for seg in segments:
        lines.append({'@x1':str(seg[0]), '@y1':str(seg[1]), '@x2':str(seg[2]), '@y2':str(seg[3]), '@stroke-width':str(width), '@stroke':color})
    return lines


def polygon_to_lines(polygon):
    # the first node in polygon is the upper-rightmost tip; append the corner preceding it so that we can process the tip correctly.
    polygon.insert(0,polygon[len(polygon)-1])
    polygon.insert(0,polygon[len(polygon)-2])

    global max_x, min_x, points
    
    (root_level, otu_level) = set_otu_and_root_level(polygon)
#     points = []
    lines = []
    lines.append([polygon[len(polygon)-1][0], polygon[len(polygon)-1][1], polygon[0][0], polygon[0][1]])
    last_node = polygon[0]
    for i in range(1, len(polygon)-1):
        node = polygon[i]
        line = [last_node[0], last_node[1], node[0], node[1]]
        # don't add zero-length lines
        if (line[0] == line[2]) and (line[1] == line[3]):
            continue
        else:
            lines.append(line)
        last_node = node

    horiz_line_set = set()
    vert_line_set = set()
    for line in lines:
        # sanity check lines:
        # if (math.fabs(line[2]-line[0]) < 2) and (math.fabs(line[3]-line[1]) < 2):
            # print "hey, weird line %s" % line
    
        forward_line = (line[0],line[1],line[2],line[3])
        if (line[0] == line[2]): # vertical lines
            if (line[1] < line[3]):
                forward_line = (line[0],line[3],line[2],line[1])
            else:
                vert_line_set.add('%d %d %d %d' % forward_line)
        elif (line[1] == line[3]): #horiz line
            if (line[0] > line[2]):
                forward_line = (line[2],line[1],line[0],line[3])
            else:
                horiz_line_set.add('%d %d %d %d' % forward_line)
    
    # sort all the horiz lines by y-value:
    horiz_lines_list = []
    for line in horiz_line_set:
        coord = line.split(' ')
        x1 = int(coord[0])
        y1 = int(coord[1])
        x2 = int(coord[2])
        y2 = int(coord[3])
        
        horiz_lines_list.append([x1, y1, x2, y2])
        
    horiz_lines_binned = sort_lines(horiz_lines_list, 1, 0)

    vert_lines_list = []
    # adjust the y2 of the vertical lines to match a horizontal line
    for line in vert_line_set:
        coord = line.split(' ')
        x = int(coord[0])
        y1 = int(coord[1])
        y2 = int(coord[3])

        my_nodeline = None
        for bin in horiz_lines_binned:
            if (y1 > bin[0][1]):
                continue
            else:
                for nodeline in bin:
                    if (nodeline[0] <= x) and (x <= nodeline[2]):
                        my_nodeline = nodeline
                        y1 = nodeline[1]
                        break
                if my_nodeline is not None:
                    break
        line = [x, y1, x, y2]

        vert_lines_list.append(line)

    # sort all the vertical lines by x-value:
    vert_lines_binned = sort_lines(vert_lines_list, 0, 1)

    # coalesce the vertical nodes
    for i in range(len(vert_lines_binned)):
        bin = vert_lines_binned[i]
        new_bin = []
        curr_node = bin[0]
        for j in range(1,len(bin)):
            line = bin[j]
            if (curr_node[1] == line[3]):
                curr_node[1] = line[1]
            else:
                new_bin.append(curr_node)
                curr_node = bin[j]
        new_bin.append(curr_node)
        vert_lines_binned[i] = new_bin
        
    # straighten out the horizontal lines:
    horiz_lines = []
    node_lines = []
    changes_made = True
    while changes_made:
        changes_made = False
        fix_bin = []
        for line in horiz_line_set:
            coord = line.split(' ')
            x1 = int(coord[0])
            y = int(coord[1])
            x2 = int(coord[2])

            # if x2 is at otu_level, we don't need to worry about that end.
            if x2 < otu_level:        
                my_nodeline = None
                i = 0
                while i < len(vert_lines_binned):
                    bin = vert_lines_binned[i]
                    if len(bin) == 0:
                        i += 1
                        continue
#                     print "%s %s" % (line, str(bin[0]))
                    if x2 > bin[0][0]:
                        i += 1
                        continue
                    else:
                        for nodeline in bin:
                            if (nodeline[3] < y) and (y < nodeline[1]):
                                my_nodeline = nodeline
                                x2 = nodeline[0]
                                break
                            elif y == nodeline[1]:
                                bin.remove(nodeline)
                                x2 = vert_lines_binned[i-1][0][0]
                                horiz_line_set.remove(line)
                                horiz_line_set.add('%d %d %d %d' % (x1, y, x2, y))
                                my_nodeline = nodeline
                                changes_made = True
                                break
                            elif (y == nodeline[3]):
                                # this is a problematic nodeline: 
                                # edges shouldn't lead into the corner of a node.
                                # print "%s removing %s" % (str([x2,y]),nodeline)
                                x2 = vert_lines_binned[i-1][0][0]
                                horiz_line_set.remove(line)
                                horiz_line_set.add('%d %d %d %d' % (x1, y, x2, y))
                                nodeline[3] = (nodeline[3] - nodeline[1]) + nodeline[3] -5
                                # print "nodeline is now %s" % str(nodeline)
                                my_nodeline = nodeline
                                changes_made = True
                                break
                        if my_nodeline is not None:
                            node_lines.append(my_nodeline)
                            break
                        i += 1
            # if x1 is at root_level, we don't need to worry about that end.
            if x1 > root_level and not changes_made:    
                my_nodeline = None
                i = 0
                while i < len(vert_lines_binned):
                    bin = vert_lines_binned[i]
                    if len(bin) == 0:
                        i += 1
                        continue
#                     print "%s %s" % (str([x1,y]), str(bin[0]))
                    if (x1 > bin[0][0]):
                        i += 1
                        continue
                    else:                    
                        for nodeline in vert_lines_binned[i]:
                            if (nodeline[3] < y) and (y < nodeline[1]):
                                my_nodeline = nodeline
                                x1 = nodeline[0]
                                break
                            elif (y == nodeline[1]):
#                                 print "%s modifying %s" % (str([x1,y]),nodeline)
                                x1 = nodeline[0]
                                my_nodeline = nodeline
                                break
                            elif y == nodeline[3]:
#                                 points.append([x1,y])
                                # print "%s touches %s" % (str([x1,y]),nodeline)
                                my_nodeline = nodeline
                                x1 = nodeline[0]
                                break
                        if my_nodeline is not None:
                            node_lines.append(my_nodeline)
                            i += 1
                            break
                    if my_nodeline is None:
                        x1 = root_level
                    i += 1
            line = [x1,y,x2,y]
            horiz_lines.append(line)
    lines = []
    for line in node_lines:
            lines.append([line[0],line[3],line[2],line[1]])
    
    lines.extend(horiz_lines)
    return lines

def sort_lines(lines, key1, key2):
    bin_by_x1 = {}
    for line in lines:
        if line[key1] not in bin_by_x1:
            bin_by_x1[line[key1]] = []
        bin_by_x1[line[key1]].append(line)
    
    bin_keys = list(bin_by_x1.keys())
    bin_keys.sort()
    
    final_lines = []
    for bin in bin_keys:
        bin_by_x1[bin].sort(cmp=lambda x,y: cmp(x[key2], y[key2]))
        final_lines.append(bin_by_x1[bin])
    return final_lines
       
def set_otu_and_root_level(polygon):
    # find the maximum x-value
    global otu_level, root_level, max_x, min_x
    root_level = max_x
    otu_level = min_x

    for node in polygon:
        if (node[0] > otu_level):
            otu_level = node[0]
        if (node[0] < root_level):
            root_level = node[0]
    return (otu_level, root_level)

# remove all in-between singletons from a polygon
def straighten_polygon(polygon):
    print "straightening polygon"
    changes_made = False
    # for convenience:
    x = 0
    y = 1
    global points
    points = []
    extend_polygon(polygon)
#     print "starting as: %s " % str(polygon)
    start_node = polygon[0]
    index = 2
    while (polygon[index] != start_node):
        # we're looking at an array of nodes: nodes[2] is equivalent to polygon[index]
        nodes = polygon[index-2:index+3]
        message = "looking at %s\n" % str(nodes)

        remove_index_node = False
        
        is_straight = False
        is_right = False
        is_tip = False
        if (nodes[1][x] == nodes[2][x]): # if first two are vertical:
            if (nodes[3][x] == nodes[2][x]):
                message += "  vertical line "
                is_straight = True
            elif (nodes[3][y] == nodes[2][y]):
                message += "  right angle\n"
                is_right = True
        elif (nodes[1][y] == nodes[2][y]): # if first two are horizontal:
            if (nodes[3][y] == nodes[2][y]):
                message += "  horiz line1 "
                is_straight = True
                if (nodes[3][x] <= nodes[2][x]) and (nodes[1][x] <= nodes[2][x]):
                    message += "  tip1 " + str(nodes[2])
                    is_tip = True
            elif (nodes[3][x] == nodes[2][x]):
                message += "  right angle\n"
                is_right = True
        elif (nodes[2][x] == nodes[3][x]): # if second two are vertical:
            if (nodes[1][x] == nodes[2][x]):
                message += "  vertical line "
                is_straight = True
            elif (nodes[1][y] == nodes[2][y]):
                message += "  right angle\n"
                is_right = True
        elif (nodes[3][y] == nodes[2][y]): # if second two are horizontal:
            if (nodes[1][y] == nodes[2][y]):
                message += "  horiz line2 "
                is_straight = True
            elif (nodes[1][x] == nodes[2][x]):
                message += "  right angle\n"
                is_right = True

        if is_right or is_tip:
            remove_index_node = False
        elif is_straight and not is_tip:
            message += "  is straight\n"
            remove_index_node = True
        else:
            res = ""
            # law of cosines:
            # we want the angle internal to node1.
            # if a, b, c are the lengths of sides of the triangle described by nodes 1,2,3
            # the angle of node1 is acos((a^2 + b^2 - c^2)/(2*a*b))
            # alternatively, if we keep the squares, acos((a2 + b2 - c2)/(2*sqrt(a2*b2)))
            # a,b,c can be solved with the Pythagorean Theorem:
            a2 = ((nodes[2][x]-nodes[1][x]) * (nodes[2][x]-nodes[1][x])) + ((nodes[2][y]-nodes[1][y]) * (nodes[2][y]-nodes[1][y]))
            b2 = ((nodes[2][x]-nodes[3][x]) * (nodes[2][x]-nodes[3][x])) + ((nodes[2][y]-nodes[3][y]) * (nodes[2][y]-nodes[3][y]))
            c2 = ((nodes[3][x]-nodes[1][x]) * (nodes[3][x]-nodes[1][x])) + ((nodes[3][y]-nodes[1][y]) * (nodes[3][y]-nodes[1][y]))
            
            theta = math.acos((a2 + b2 - c2)/(2*math.sqrt(a2*b2)))

            # is their angle really flat
            if 180 - math.degrees(theta) < 30:
                remove_index_node = True
                res = "basically straight"
            
            message += "  theta is %f: %s\n" % (math.degrees(theta), res)

        if remove_index_node:
            # remove node 2 (which is polygon[index]), don't increment:
            print message + "  removing node %s" % str(nodes[2])
            polygon.pop(index)
            index = index - 1
            changes_made = True

        # end loop by incrementing index
        index = index + 1
    
    trim_polygon(polygon)
    if changes_made:
        straighten_polygon(polygon)
#     print "ending as: %s " % str(polygon)
    return


def pointify_tips(polygon):
    print "pointifying tips of polygon"
    extend_polygon(polygon)
    changes_made = False
    # for convenience:
    x = 0
    y = 1
    global points
    points = []
    start_node = polygon[0]
    index = 2
    while (polygon[index] != start_node):
        # we're looking at an array of nodes: nodes[2] is equivalent to polygon[index]
        nodes = polygon[index-2:index+3]
        
        remove_index_node = False

        # tips can either already be pointy or be blunt.
        
        # if nodes[0] and nodes[1] increases in x, we're looking for leaf tips.
        if nodes[1][x] >= nodes[0][x]:
            print "  examining %s" % str(nodes)
            # a tip has to have the y vals of nodes[2] be greater than that of nodes[1]
            if nodes[2][y] <= nodes[1][y]:
                print "    looking for tip in %s" % str(nodes)
                # if nodes[2] goes back (has smaller x than nodes[1]), it's a pointy tip. Leave it alone.
                if nodes[2][x] < nodes[1][x]:
                    print "    pointy tip at %s" % str(nodes[1])
                    remove_index_node = False
                # if nodes[2] still increases in x but nodes[3] is smaller, it's a blunt tip: remove nodes[2], the index node.
                elif nodes[3][x] <= nodes[2][x] <= nodes[1][x]:
                    remove_index_node = True

        if remove_index_node:
            # remove node 1 (which is polygon[index-1]), don't increment:
            print "    looking at %s" % str(nodes)
            print "      removing node %s" % str(nodes[1])
            polygon.pop(index-1)
            index = index - 1
            changes_made = True

        # end loop by incrementing index
        index = index + 1
    
    trim_polygon(polygon)
#     print "ending as: %s " % str(polygon)
    return

# rotate the polygon so that it starts at the smallest x,y node
def rotate_polygon(polygon):
    smallest_node = polygon[0]
    smallest_index = 0
    i = 1
    while (i < len(polygon)):
        if polygon[i] < smallest_node:
            smallest_node = polygon[i]
            smallest_index = i
        i = i + 1
    print "smallest node is %s" % str(smallest_node)
        
    # add front half to the end of the list
    polygon.extend(polygon[:smallest_index])
    # delete the first half of the list
    del polygon[:smallest_index]
        
def remove_duplicate_points(polygon):
    print "removing dups in polygon: starts as %s, ends as %s" %(polygon[:3],polygon[len(polygon)-10:])
    working_set = list(polygon)
    polygon[:] = [working_set.pop(0)]
    while len(working_set) > 2:
        node1 = polygon[len(polygon) - 1]
        node2 = working_set.pop(0)
        if node1 != node2:
            polygon.append(node2)
            working_set.insert(0, node2)
    print "done removing dups in polygon: starts as %s, ends as %s" %(polygon[:3],polygon[len(polygon)-10:])


def straighten_horizontal_lines(polygon):
    print "straightening horizontal lines"
    extend_polygon(polygon)
    index = 2
    start_node = polygon[0]
    y = 1
    while (polygon[index] != start_node):
        # we're looking at an array of nodes: nodes[2] is equivalent to polygon[index]
        nodes = polygon[index-2:index+3]

        # look at the y-values of these: are they on a horizontal line?
        average_y = (nodes[0][y] + nodes[1][y] + nodes[2][y] + nodes[3][y] + nodes[4][y]) / 5
        print "  straightening along %s: looking at %s" % (average_y, str(nodes))
        horizontal = True
        for i in range(5):
            if (nodes[i][y] - average_y) > 5:
                horizontal = False

        if horizontal:
            print "  straight horizontal line along %s" % str(average_y)
            for i in range(5):
                nodes[i][y] = average_y
            print "  straight horizontal line along %s, straightened:  %s" % (average_y, str(nodes))
        # end loop by incrementing index
        index = index + 1
    trim_polygon(polygon)


def trim_polygon(polygon):
    remove_duplicate_points(polygon)
    print "trimming polygon: starts as %s, ends as %s" %(polygon[:3],polygon[len(polygon)-10:])
    # trim any extra nodes off the end:
    while polygon[0] != polygon[len(polygon)-1]:
        lastnode = polygon.pop()  
        print "  trimmed %s" % (str(lastnode))
    # #pop one more:
    lastnode = polygon.pop()
    print "finally trimmed %s" % (str(lastnode))
    rotate_polygon(polygon)


def extend_polygon(polygon):
    rotate_polygon(polygon)
    polygon.extend(polygon[0:4])
    print "extending polygon: starts as %s, ends as %s" %(polygon[:3],polygon[len(polygon)-10:])


def find_otus(polygon):
    # for convenience:
    x = 0
    y = 1
#     print polygon

    # the first node in the polygon is the upper-rightmost tip
    global points
    points = []
    polygon.insert(0,polygon[len(polygon)-1])
    
    # we need to make sure we start with the last thing in polygon
    new_polygon = []
    new_polygon.append(polygon.pop())
    new_polygon.append(polygon.pop(0))
    new_polygon.append(polygon.pop(0))

    while len(polygon) >= 0:
        node2 = new_polygon.pop()
        node1 = new_polygon.pop()
        node0 = new_polygon.pop()

        if (node1[x] > node2[x]) and (node1[x] > node0[x]):
            points.append(node1)
#             print node1

        #### FINALLY: append nodes
        new_polygon.append(node0)
        new_polygon.append(node1)
        new_polygon.append(node2)
        if len(polygon) == 0:
            break
        
        new_polygon.append(polygon.pop(0))
    new_polygon.pop(0)
#     print new_polygon
    return new_polygon

def path_to_polygon(path):
    polygon = []
    global max_x, max_y, min_x, min_y
    new_path = Path()    
    for segment in parse_path(path):
        new_path.append(Line(segment.start, segment.end))
    new_path.closed = True
    raw_path = new_path.d()
    
    # search for compound path bits and remove them
    path_bits = re.findall('M.+?[ZM]', raw_path)
    if len(path_bits) > 0:
        raw_path = path_bits[0]
        raw_path_size = len(re.findall(',',raw_path))
        for bit in path_bits:
            bit_size = len(re.findall(',',bit))
            if bit_size > raw_path_size:
                raw_path = bit
                raw_path_size = bit_size
    
    # convert to simple list of nodes
    nodes = re.findall('[ML]\s*(\d+\.*\d*,\d+\.*\d*)\s*', raw_path)
    for n in nodes:
        coords = n.split(',')
        if max_x < int(coords[0]):
            max_x = int(coords[0])
        if max_y < int(coords[1]):
            max_y = int(coords[1])
        if min_x > int(coords[0]):
            min_x = int(coords[0])
        if min_y > int(coords[1]):
            min_y = int(coords[1])
        polygon.append([int(coords[0]), int(coords[1])])
    polygon.pop()
    return polygon
    
def nodes_to_path(nodes):
    path_points = []
    for point in nodes:
        path_points.append('%d %d' % (point[0],point[1]))
    return 'M' + 'L'.join(path_points) + 'Z'
          
def nodes_to_circles(nodes):
    circlelist = []
    for i in range(len(nodes)):
        circledict = {}
        coords = nodes[i]
        circledict['@r'] = '3'
        circledict['@stroke'] = 'black'
        circledict['@stroke-width'] = '1'
        circledict['@fill'] = 'yellow'
        circledict['@cx'] = str(coords[0])
        circledict['@cy'] = str(coords[1])
        circlelist.append(circledict)
    return circlelist

def threshold_area(rect, threshold):
    mins = rect[0]
    maxs = rect[2]
    min_x = float(mins[0])
    min_y = float(mins[1])
    max_x = float(maxs[0])
    max_y = float(maxs[1])
    if abs(float(max_y - min_y) / float(total_height)) > float(threshold):
        return True
    if abs(float(max_x - min_x) / float(total_width)) > float(threshold):
        return True
    return False
    
def bounding_box(polygon):
    x_points = []
    y_points = []
    for point in polygon:
        coord = point
        x_points.append(coord[0])
        y_points.append(coord[1])
    max_x = max(x_points)
    max_y = max(y_points)
    min_x = min(x_points)
    min_y = min(y_points)
    return [[min_x, min_y],[min_x, max_y],[max_x, max_y],[max_x, min_y]]

def define_text_boxes(boxpaths, nodes, segments, otus):
    edges = []
    for seg in segments:
        [x1, y1, x2, y2] = seg
        if (y1 == y2):
            edges.append(seg)
            
    otu_edges = []
    branch_edges = []
    otu_label_dict = {}
    branch_label_dict = {}
    unknown_label_dict = {'x':[]}

    # find edges that belong to otus:
    for edge in edges:
        [x1, y1, x2, y2] = edge
        if [x2, y1] in otus:
            otu_edges.append(edge)
            otu_label_dict['%d %d %d %d' % (edge[0],edge[1],edge[2],edge[3])] = []
        else:
            branch_edges.append(edge)
            branch_label_dict['%d %d %d %d' % (edge[0],edge[1],edge[2],edge[3])] = []

    for box in boxpaths:
        [box_x1, box_y1] = box[0]
        [box_x2, box_y2] = box[2]
        marked = False
        # if text occurs above or below a branch_edge, it's a branch label
        for edge in branch_edges:
            [edge_x1, edge_y1, edge_x2, edge_y2] = edge
            # sometimes labels are longer than the branch:
            edge_x1 -= 10
            edge_x2 += 10
            if (box_x1 >= edge_x1) and (box_x2 <= edge_x2):
                # above the branch:
                # extend the y-boundaries of the edge a little:
                edge_y1 = edge[1] - 15
                edge_y2 = edge[1] - 5
                if (box_y1 <= edge_y2) and (box_y2 >= edge_y1):
                    branch_label_dict['%d %d %d %d' % (edge[0],edge[1],edge[2],edge[3])].extend(box)
                    marked = True
                # below the branch:
                # extend the y-boundaries of the edge a little:
                edge_y1 = edge[1] + 5
                edge_y2 = edge[1] + 15
                if (box_y1 <= edge_y2) and (box_y2 >= edge_y1):
                    branch_label_dict['%d %d %d %d' % (edge[0],edge[1],edge[2],edge[3])].extend(box)
                    marked = True
        
        if marked == False:
            # if text occurs to the right of an otu_edge's x2, it's an otu_label
            for edge in otu_edges:
                [edge_x1, edge_y1, edge_x2, edge_y2] = edge
                # extend the y-boundaries of the edge a little:
                edge_y1 -= 2
                edge_y2 += 2
                if (box_y1 <= edge_y2) and (box_y2 >= edge_y1) and (box_x1 > edge_x2):
                    otu_label_dict['%d %d %d %d' % (edge[0],edge[1],edge[2],edge[3])].extend(box)
                    marked = True
        
        if marked == False:
            unknown_label_dict['x'].append(box)
    
    for edge in otu_label_dict.keys():
        if len(otu_label_dict[edge]) > 0:
            [x1,y1,x2,y2] = re.split(' ', edge)
            large_box = bounding_box(otu_label_dict[edge])
            otu_label_dict[str([int(x2),int(y2)])] = large_box
        otu_label_dict.pop(edge, None)

    for edge in branch_label_dict.keys():
        if len(branch_label_dict[edge]) > 0:
            [x1,y1,x2,y2] = re.split(' ', edge)
            large_box = bounding_box(branch_label_dict[edge])
            branch_label_dict[str([int(x2),int(y2)])] = large_box
        branch_label_dict.pop(edge, None)

    return (otu_label_dict, branch_label_dict, unknown_label_dict)

if __name__ == '__main__':
    main()
