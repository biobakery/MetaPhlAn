#from dendropy import Tree
from dendropy.datamodel.treemodel import Tree
from math import sqrt
try:
    from queue import Queue # python 3
except:
    from Queue import Queue # python 2
from copy import deepcopy
from sys import stdout
from .Tree_extend import Centroid_Tree, MPR_Tree


class TreeInduced:
    def __init__(self,bestLCA=None):
        self.bestLCA = bestLCA
        self.records = {}

class Entry:
    def __init__(self,bestLCA=None):
        self.level = 0
        self.info = TreeInduced(bestLCA=bestLCA)
        self.backtrack = None
        self.removed = None
        self.retained = None

class TreeFilter:       
    def __init__(self, ddpTree = None, tree_file = None, scaling = None,  schema = "newick", centroid_reroot = False):
        a_tree = Centroid_Tree(ddpTree=ddpTree,tree_file=tree_file,schema=schema)

        if centroid_reroot:
#            print("Rerooting at centroid ...")
            a_tree.Reroot()
#        if tree_file:
#            self.ddpTree = Tree.get_from_path(tree_file,schema,preserve_underscores=True)
#        else:  
#            self.ddpTree = ddpTree

        self.ddpTree = a_tree.ddpTree

        self.bestLCA = None
        self.nleaf = 0
        self.records = {}
        self.a = scaling[0]
        self.b = scaling[1]
        #print("Using a= %d , b= %d" %(self.a, self.b))

        diam = -1

        for node in self.ddpTree.postorder_node_iter():
               if node.is_leaf():
                   self.nleaf += 1
               self.__updateNode__(node)
               if (self.records[node][3] and self.records[node][3] > diam):
                   diam = self.records[node][3]
                   self.bestLCA = node

        #print(self.records[self.bestLCA][0].taxon.label)
        #print(self.records[self.bestLCA][1].taxon.label)

        #self.myQueue = [first_entry]
        self.myQueue = Queue()
        self.best_entries = []
        self.min_diams = []
    
    def __updateNode__(self,node,records=None):
        if records is None:
            records = self.records

        if node.is_leaf():
            records[node] = [node,node,0,0]
            return

        records[node] = None

        max1 = -1
        max2 = -1
        anchor1 = None
        anchor2 = None

        for ch in node.child_node_iter():
               ch_record = records[ch] if ch in records else self.records[ch]
               if ch_record is None:
                    # this child node was removed!
                    continue
               l = ch_record[2] + ch.edge_length
               if l > max1:
                   max2 = max1
                   max1 = l
                   anchor2 = anchor1
                   anchor1 = ch_record[0]
               elif l > max2:
                   max2 = l
                   anchor2 = ch_record[0]

        MAX = max1+max2 if anchor2 else None
        records[node] = [anchor1,anchor2,max1,MAX] if anchor1 else None

    def __substitute_anchor__(self,entry,anchor,inherit_info=True):
        #print("previous level: " + str(entry.level))
        #print("previous anchors: " + self.__get_anchor1__(entry).taxon.label,self.__get_anchor2__(entry).taxon.label)
        anchor1 = self.__get_anchor1__(entry)
        anchor2 = self.__get_anchor2__(entry)

        if anchor is anchor1:
            retained_anchor = anchor2
        elif anchor is anchor2:
            retained_anchor = anchor1
        else:
            print("ERROR : TreeFilter.__substitute_anchor__: The request anchor to be removed is not one of the two anchors!")
            return None
       
        new_entry = Entry()
        if inherit_info:
            new_entry.info = entry.info
        else:
            #new_entry.info = deepcopy(entry.info)
            new_entry.info.bestLCA = entry.info.bestLCA
            new_entry.info.records = {}
            for key in entry.info.records:
                new_entry.info.records[key] = entry.info.records[key]
        
        new_entry.level = entry.level + 1

        # update node records from anchor to root
        new_entry.info.records[anchor] = None 
        curr_node = anchor.parent_node
        new_best_LCA = None
        new_best_record = None

        while curr_node:
            self.__updateNode__(curr_node,records=new_entry.info.records)
            node_record = self.__lookup__(new_entry,curr_node)
            if node_record is not None and node_record[1] is not None and node_record[3] is not None and (new_best_LCA is None or node_record[3] > new_best_record[3]): 
                new_best_LCA = curr_node
                new_best_record = node_record
            curr_node = curr_node.parent_node

        # lookup for the new bestLCA: it must be on the path from the retained anchor to the current bestLCA of the entry
        curr_node = retained_anchor.parent_node
        node_record = self.__lookup__(new_entry,curr_node)

        while curr_node is not entry.info.bestLCA:
            if node_record is not None and node_record[1] is not None and node_record[3] is not None and (new_best_LCA is None or node_record[3] > new_best_record[3]):
                new_best_LCA = curr_node
                new_best_record = node_record
            curr_node = curr_node.parent_node        
            node_record = self.__lookup__(new_entry,curr_node)
         
        if node_record is not None and node_record[1] is not None and node_record[3] is not None and (new_best_LCA is None or node_record[3] > new_best_record[3]):
            new_best_LCA = curr_node

        new_entry.info.bestLCA = new_best_LCA

        # backtrack
        new_entry.backtrack = entry
        new_entry.removed = anchor
        new_entry.retained = retained_anchor

        #print(new_entry.level)

        #print(new_entry.removed.taxon.label + " removed")
        #print("current anchors: " + self.__get_anchor1__(new_entry).taxon.label,self.__get_anchor2__(new_entry).taxon.label)

        return new_entry

    def __substitute_anchor1__(self,entry,inherit_info=True):
        return self.__substitute_anchor__(entry,self.__get_anchor1__(entry),inherit_info=inherit_info)

    def __substitute_anchor2__(self,entry,inherit_info=True):
        return self.__substitute_anchor__(entry,self.__get_anchor2__(entry),inherit_info=inherit_info)
   
    def __get_diam__(self,entry):
         try:
            d = entry.info.records[entry.info.bestLCA][3]
         except:
            d = self.records[entry.info.bestLCA][3]
         return d

    def __get_anchors__(self,entry):
        if entry.info.bestLCA in entry.info.records:
            anchor1 = entry.info.records[entry.info.bestLCA][0]
            anchor2 = entry.info.records[entry.info.bestLCA][1]
        else:
            anchor1 = self.records[entry.info.bestLCA][0]
            anchor2 = self.records[entry.info.bestLCA][1]
    
        if anchor1 is entry.retained:
            temp = anchor2 
            anchor2 = anchor1
            anchor1 = temp

        return anchor1, anchor2

    def __get_anchor1__(self,entry):
        if entry.info.bestLCA in entry.info.records:
            node = entry.info.records[entry.info.bestLCA][0]
        else:
            node = self.records[entry.info.bestLCA][0]
        return node


    def __get_anchor2__(self,entry):
        if entry.info.bestLCA in entry.info.records:
            node = entry.info.records[entry.info.bestLCA][1]
        else:
            node = self.records[entry.info.bestLCA][1]
        return node
    
    def __lookup__(self,entry,node):
        if node in entry.info.records:
            node_record = entry.info.records[node]
        else:
            node_record = self.records[node]

        return node_record
    

    def __default_d__(self,DEFAULT_MIN=0):
        return min(self.nleaf//self.a,max(DEFAULT_MIN,int(self.b*sqrt(self.nleaf))))

    def optFilter(self,d=None):
        d = self.__default_d__() if d is None else d

        print("Solving k-shrink with k = " + str(d));

        first_entry = Entry(bestLCA=self.bestLCA)
        #print(self.__get_anchor1__(first_entry).taxon.label,self.__get_anchor2__(first_entry).taxon.label)
        self.myQueue.put(first_entry)
        
        curr_level = -1
        min_diam = None 
        best_entry = None

        while 1:
            if self.myQueue.empty():
                break
            curr_entry = self.myQueue.get()
            if curr_entry.info.bestLCA is None:
                break
            diam = self.__get_diam__(curr_entry)
            
            anchor1,anchor2 = self.__get_anchors__(curr_entry)
            if curr_entry.level != curr_level:
                # reached a next level
                curr_level = curr_entry.level
                if min_diam is not None:
                    self.min_diams.append(min_diam)
                    self.best_entries.append(best_entry)
                    
                    #if best_entry:
                        #print("Best: ")    
                        #print(self.__lookup__(best_entry,best_entry.info.bestLCA)[0].taxon.label)
                        #print(self.__lookup__(best_entry,best_entry.info.bestLCA)[1].taxon.label)
                        #print(self.__lookup__(best_entry,best_entry.info.bestLCA)[3])
                    #print(curr_entry.level)    
                    #print(self.__lookup__(curr_entry,curr_entry.info.bestLCA)[0].taxon.label)
                    #print(self.__lookup__(curr_entry,curr_entry.info.bestLCA)[1].taxon.label)
                    #print(self.__lookup__(curr_entry,curr_entry.info.bestLCA)[3])
                min_diam = diam
                best_entry = curr_entry
                if curr_level < d:
                    # add 2 entries
                    self.myQueue.put(self.__substitute_anchor__(curr_entry,anchor1,inherit_info=False))
                    self.myQueue.put(self.__substitute_anchor__(curr_entry,anchor2,inherit_info=True))
            else:
                #print(curr_entry.level)
                #print(self.__lookup__(curr_entry,curr_entry.info.bestLCA)[0].taxon.label)
                #print(self.__lookup__(curr_entry,curr_entry.info.bestLCA)[1].taxon.label)
                #print(self.__lookup__(curr_entry,curr_entry.info.bestLCA)[3])
                if diam < min_diam:
                    min_diam  = diam
                    best_entry = curr_entry
                if curr_level < d:
                    # add 1 entry
                    #if anchor1 is curr_entry.retained:
                    #    anchor = anchor2
                    #elif anchor2 is curr_entry.retained:
                    #    anchor = anchor1
                    #else:
                    #    print("FilterTree.optFilter: anchors inconsistent between child and parent!")
                    #    return False
                    self.myQueue.put(self.__substitute_anchor__(curr_entry,anchor1,inherit_info=True))

            curr_entry.info = None
        self.min_diams.append(min_diam)
        self.best_entries.append(best_entry)


        #print(self.min_diams)

    def __prune_taxon__(self,taxon):
        pnode = taxon.parent_node
        pnode.remove_child(taxon)
        if pnode.num_child_nodes() < 2:
            if pnode is self.ddpTree.seed_node:
            # remove pnode and make its only child the new root
                ch = pnode.child_nodes()[0]
                pnode.remove_child(ch)
                self.ddpTree.seed_node = ch
            else:
            # remove pnode and make its only child the child of its parent
                gnode = pnode.parent_node
                gnode.remove_child(pnode)
                gnode.add_child(pnode.child_nodes()[0])

    def filterOut(self,d=None,fout=None):
        d = self.__default_d__() if d is None else d
        entry = self.best_entries[d]
        while entry.backtrack is not None:
            if fout:
                fout.write(entry.removed.taxon.label + "\t")
            self.__prune_taxon__(entry.removed)
            entry = entry.backtrack

        return self.ddpTree

    def list_removals_reverse(self,d=None,fout=None): 
        d = self.__default_d__() if d is None else d
        last_entry = self.best_entries[d]
        rm_list = []
        def __list__(entry):
            if entry.backtrack is not None:
                __list__(entry.backtrack)
                if fout:
                    fout.write(entry.removed.taxon.label + " ")
                rm_list.append(entry.removed.taxon.label)

        __list__(last_entry)
        return rm_list

    def list_removals(self,d=None,fout=None):
        d = self.__default_d__() if d is None else d
        entry = self.best_entries[d]
        rm_list = []
        while entry.backtrack is not None:
            if fout:
                fout.write(entry.removed.taxon.label + " ")
            rm_list.append(entry.removed.taxon.label)
            entry = entry.backtrack

        return rm_list
