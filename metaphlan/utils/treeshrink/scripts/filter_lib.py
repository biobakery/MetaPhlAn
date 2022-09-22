import operator
from dendropy import Tree
from copy import deepcopy
from .Tree_extend import MPR_Tree, MV00_Tree,MVDF_Tree

def filter_branch(a_tree,root_method=None,unit_length=None,low_percentile=0,high_percentile=1,factor=1):
    a_tree.deroot()
    branch_list = list_branch(a_tree)
    d = estimate_diameter(a_tree,branch_list,root_method=root_method,unit_length=unit_length,low_percentile=low_percentile,high_percentile=high_percentile)
    thres = d*factor
    print("Branch length threshod: ",thres)
    count_leaves(a_tree)
    for br in branch_list:
        #print(br.length)
        if br.length > thres:
            remove_branch(a_tree,br)

def count_leaves(a_tree):
    for node in a_tree.postorder_node_iter():
        if node.is_leaf():
            node.nleaf = 0
        else:
            node.nleaf = sum([ch.nleaf for ch in node.child_node_iter()])


def remove_branch(a_tree,br):
            p = br.tail_node
            c = br.head_node
            if c.nleaf > a_tree.seed_node.nleaf/2:
                p.remove_child(c)
                a_tree.seed_node = c
            elif p is not None:
                p.remove_child(c,suppress_unifurcations=True)
                while p and p.num_child_nodes() == 0:
                    c = p
                    p = p.parent_node
                    if p is not None:
                        p.remove_child(c,supress_unifurcations=True)


def list_branch(a_tree):
    branch_list = []
    for br in a_tree.preorder_edge_iter():
        if br.tail_node is not None:
            branch_list.append(br)
    branch_list.sort(key=lambda x: x.length)
    
    return branch_list


def estimate_diameter(a_tree,branch_list,root_method=None,unit_length=None,low_percentile=0,high_percentile=1):
    # root_method: 'MV00','MVDF', or 'MP'
    # unit_length: use the specified length to be the branch length of a unit-branch tree 
    # special options: 'median': use the median branch length as the length for the unit tree
    #                  'avg': use the average instead of the median
    # low_percentile: the lowest accepted percentile of branches (branches shorter than this will be set to this value)
    # high_percentile: the highest accepted percentile of branches ((branches longer than this will be set to this value)

    def __brlen(br,unit_based,low_thres,high_thres):
        if unit_based:
            return 1
        elif br < low_thres:
            return low_thres
        elif br > high_thres:
            return high_thres
        else:
            return br
    
    def __compute_max_distance(unit_based=False,low_thres=-1,high_thres=10**8):
        max_distance = 0

        for node in a_tree.postorder_node_iter():
            if node.is_leaf():
                node.max_br_below = 0
            else:
                children = node.child_nodes()
                l1 = children[0].max_br_below + __brlen(children[0].edge_length,unit_based,low_thres,high_thres) 
                l2 = children[1].max_br_below + __brlen(children[1].edge_length,unit_based,low_thres,high_thres) 
                max1 = max(l1,l2)
                max2 = min(l1,l2)                
                
                i = 2
                while i < len(children):
                    l = children[i].max_br_below + __brlen(children[i].edge_length,unit_based,low_thres,high_thres) 
                    if l > max1:
                        max1 = l 
                    elif l > max2:
                        max2 = l
                    i += 1

                max_distance = max(max_distance,max1 + max2)
                node.max_br_below = max1
        
        return max_distance


    def __unit_based_diameter(unit='median'):
        if unit=='median':    
            unit_length = branch_list[int(len(branch_list)*0.5)].length
        elif unit=='avg':
            unit_length = 0
            for br in branch_list:
                unit_length += br.length
                unit_length /= len(branch_list)
        else:
            unit_length = unit

        return unit_length* __compute_max_distance(unit_based=True)

    def __percentile_based_diameter(low_percentile=0,high_percentile=1):
        low_thres = branch_list[int(len(branch_list)*low_percentile)].length
        high_thres = branch_list[int(len(branch_list)*high_percentile)].length
        return __compute_max_distance(low_thres=low_thres,high_thres=high_thres)

    if root_method:
        if root_method == 'MP':
            rooted_tree = MPR_Tree(ddpTree=a_tree)
        elif root_method == 'MV00':
            rooted_tree = MV00_Tree(ddpTree=a_tree)
        else:
            rooted_tree = MVDF_Tree(ddpTree=a_tree)

        rooted_tree.Reroot()
        D = rooted_tree.compute_ingroup_distances()
        D.sort()
        print(D)
        return D[int(len(D)*0.5)]

    if unit_length is not None:
        return __unit_based_diameter(unit=unit_length)
    
    return __percentile_based_diameter(low_percentile=low_percentile, high_percentile=high_percentile)    
