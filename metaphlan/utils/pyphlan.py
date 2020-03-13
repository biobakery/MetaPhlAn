from Bio import Phylo
from Bio.Phylo import PhyloXML
from Bio.Phylo import PhyloXMLIO
from collections import defaultdict as ddict
from Bio.Phylo.PhyloXML import Property as Prop
from Bio.Phylo.PhyloXML import Clade as PClade
from Bio.Phylo.BaseTree import Tree as BTree
from Bio.Phylo.BaseTree import Clade as BClade
import string
from numpy import pi as rpi
rpi2 = 2.0*rpi
import numpy as np
import array as arr
import collections as colls
import sys
#core_test = lambda ok,tot,pr: 1.0-st.binom.sf(ok,tot,pr)

#lev_sep = "."

uncl = "?"

# Here are three functions that I'd love to see in Biopython but they
# are not there (yet).
def partial_branch_length(clade, selective_targets):
    def _partial_branch_length_( clade, selective_targets ):
        if clade.is_terminal() and clade.name in selective_targets:
            return [clade.branch_length] if clade.branch_length else [0.0]
        if not any([c.name in selective_targets for c in clade.get_terminals()]):
            return [0.0]
        ret = [0.0]
        for c in clade.clades:
            ret += [partial_branch_length( c, selective_targets)]
        ret += [clade.branch_length] if clade.branch_length else [0.0]
        return ret
    return sum( _partial_branch_length_( clade,selective_targets  )  )

def reroot(tree, new_root):
    outgroup = new_root
    outgroup_path = tree.get_path(outgroup)
    if len(outgroup_path) == 0:
        # Outgroup is the current root -- no change
        return
    prev_blen = outgroup.branch_length
    if outgroup.is_terminal():
        # Create a new root with a 0-length branch to the outgroup
        outgroup.branch_length = 0.0
        new_root = tree.root.__class__(
                branch_length=tree.root.branch_length, clades=[outgroup])
        # The first branch reversal (see the upcoming loop) is modified
        if len(outgroup_path) == 1:
            # Trivial tree like '(A,B);
            new_parent = new_root
        else:
            parent = outgroup_path.pop(-2)
            parent.clades.pop(parent.clades.index(outgroup))
            prev_blen, parent.branch_length = parent.branch_length, prev_blen
            new_root.clades.insert(0, parent)
            new_parent = parent
    else:
        # Use the given outgroup node as the new (trifurcating) root
        new_root = outgroup
        new_root.branch_length = tree.root.branch_length
        new_parent = new_root

    # Tracing the outgroup lineage backwards, reattach the subclades under a
    # new root clade. Reverse the branches directly above the outgroup in
    # the tree, but keep the descendants of those clades as they are.
    for parent in outgroup_path[-2::-1]:
        parent.clades.pop(parent.clades.index(new_parent))
        prev_blen, parent.branch_length = parent.branch_length, prev_blen
        new_parent.clades.insert(0, parent)
        new_parent = parent

    # Finally, handle the original root according to number of descendents
    old_root = tree.root
    if outgroup in old_root.clades:
        assert len(outgroup_path) == 1
        old_root.clades.pop(old_root.clades.index(outgroup))
    else:
        old_root.clades.pop(old_root.clades.index(new_parent))
    if len(old_root) == 1:
        # Delete the old bifurcating root & add branch lengths
        ingroup = old_root.clades[0]
        if ingroup.branch_length:
            ingroup.branch_length += prev_blen
        else:
            ingroup.branch_length = prev_blen
        new_parent.clades.insert(0, ingroup)
        # ENH: If annotations are attached to old_root, do... something.
    else:
        # Keep the old trifurcating/polytomous root as an internal node
        old_root.branch_length = prev_blen
        new_parent.clades.insert(0, old_root)

    tree.root = new_root
    tree.rooted = True
    return

def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2] if len(node_path) > 1 else None

def reroot_mid_fat_edge( tree, node ):
    if tree.root == node:
        return
    fat = get_parent( tree, node )
    bl = node.branch_length
    node.branch_length = bl*0.5
    new_clade = PClade(branch_length=bl*0.5, clades = [node])
    if fat:
        fat.clades = [c for c in fat.clades if c != node] + [new_clade]
        reroot( tree, new_clade )
    else:
        tree.root.clades = [new_clade] + [c for c in tree.root.clades if c != node]
        reroot( tree, new_clade)

def clades2terms( tree, startswith = None ):
    c2t = {}
    def clades2terms_rec( c ):
        if startswith:
            if c.name and c.name.startswith( startswith ):
                c2t[c] = c.get_terminals()
        else:
            c2t[c] = c.get_terminals()
        for cc in c.clades:
            clades2terms_rec(cc)
    clades2terms_rec( tree.root )
    return c2t

def dist_matrix( tree ):
    terminals = list(tree.get_terminals())
    term_names = [t.name for t in terminals]
    # can be made faster with recursion
    for n in tree.get_nonterminals():
        n.ids = set( [nn.name for nn in n.get_terminals()]  )

    dists = dict([(n,dict([(nn,0.0) for nn in term_names])) for n in term_names])

    def dist_matrix_rec( clade ):
        bl = clade.branch_length
        if clade.is_terminal():
            for t in term_names:
                if t!=clade.name:
                    dists[clade.name][t] += bl
                    dists[t][clade.name] += bl
            return
        for t1 in clade.ids:
            for t2 in terminals:
                if t2.name not in clade.ids:
                    dists[t1][t2.name] += bl
                    dists[t2.name][t1] += bl
        for c in clade.clades:
            dist_matrix_rec( c )

    dist_matrix_rec( tree.root )
    return dists



class PpaTree:

    def __load_tree_txt__( self, fn ):
        tree = Phylo.BaseTree.Tree()
        try:
            rows = [l.decode('utf-8').rstrip().split("\t")[0] for l in
                        open(fn, 'rb')]
        except IOError:
            raise IOError()

        clades = [r.split(self.lev_sep) for r in rows]

        tree = BTree()
        tree.root = BClade()

        def add_clade_rec( father, txt_tree ):
            fl = set([t[0] for t in txt_tree])
            father.clades = []
            for c in fl:
                nclade = BClade( branch_length = 1.0,
                                 name = c )
                father.clades.append( nclade )
                children = [t[1:] for t in txt_tree if len(t)>1 and t[0] == c]
                if children:
                    add_clade_rec( nclade, children )

        add_clade_rec( tree.root, clades )
        self.ignore_branch_len = 1
        return tree.as_phyloxml()


    def __read_tree__( self, fn ):
        for ff in ['phyloxml','newick','nexus',"txt"]:
            try:
                if ff in ['txt']:
                    tree = self.__load_tree_txt__( fn )
                else:
                    tree = Phylo.read(fn, ff)
                    if len(tree.root.get_terminals()) == 1:
                        raise ValueError
            except ValueError:
                continue
            except IOError:
                sys.stderr.write("Error: No tree file found: "+fn+"\n")
                raise IOError
            except Exception:
                continue
            else:
                return tree.as_phyloxml()
        sys.stderr.write("Error: unrecognized input format "+fn+"\n")
        raise ValueError


    def __init__( self, filename, warnings = False, lev_sep = "." ):
        self.lev_sep = lev_sep
        self.warnings = warnings
        if filename is None:
            self.tree = None
            return
        try:
            self.tree = self.__read_tree__(filename)
            self.add_full_paths()
        except:
            sys.exit(0)


    def core_test( self, ok, tot, pr ):
        # scipy included here for non-compatibility with scons
        import scipy.stats as st
        if pr in self.ctc and tot in self.ctc[pr] and ok in self.ctc[pr][tot]:
            return self.ctc[pr][tot][ok]
        ret = 1.0-st.binom.sf(ok,tot,pr)
        if not pr in self.ctc: self.ctc[pr] = {}
        if not tot in self.ctc[pr]: self.ctc[pr][tot] = {}
        if not ok in self.ctc[pr][tot]: self.ctc[pr][tot][ok] = ret
        return ret

    def is_core( self, clade, targs, er = 0.95 ):
        intersection = clade.imgids & targs

        len_intersection = len(intersection)

        if len(clade.imgids) >= 2 and len_intersection < 2:
           return False, 0.0, None

        add = 0
        for subclade in clade.clades:
            if uncl in subclade.name:
                out = subclade.imgids - intersection # targs
                add += len(out)
        if add and len_intersection >= add:
            len_intersection += int(round(add/1.99))

        core = self.core_test( len_intersection, clade.nterminals, er )
        if core < 0.05 or len_intersection == 0:
            return False, core, None
        nsubclades, nsubclades_absent = 0, 0
        for subclade in set(clade.get_nonterminals()) - set([clade]):
            if uncl in subclade.full_name: # full??/
                continue
            if subclade.nterminals == 1:
                nsubclades += 1 # !!!
                if len(subclade.imgids & targs) == 0:
                    nsubclades_absent += 1
                continue

            sc_intersection = subclade.imgids & targs
            sc_len_intersection = len(sc_intersection)

            sc_add = 0
            for sc_subclade in subclade.clades:
                if uncl in sc_subclade.name:
                    sc_out = sc_subclade.imgids - sc_intersection
                    sc_add += len(sc_out)
            if add and sc_len_intersection >= sc_add:
                sc_len_intersection += int(round(sc_add/1.99))

            subcore = self.core_test( sc_len_intersection, subclade.nterminals, er )
            if subcore < 0.05:
                return False, core, None
        if nsubclades > 0 and nsubclades == nsubclades_absent:
            return False, core, None
        return True, core, intersection

    def _find_core( self, terminals, er = 0.95, root_name = None, skip_qm = True ):
        #terminals_s = set(terminals)
        def _find_core_rec( clade ):
            """
            if root_name:
                #clname = lev_sep.join( [root_name]+clade.full_name.split(lev_sep)[1:] )
                #clname = lev_sep.join( clade.full_name[1:] )
                clname = clade.full_name
            else:
                clname = clade.full_name
            """
            clname = clade.full_name
            if clade.is_terminal():
                if clade.imgid in terminals:
                    #n = terminals[clade.imgid]
                    return [(clname,1,1,
                                #n,n,n,
                                1.0)]
                return []
            if skip_qm and clade.name and uncl in clade.name:
                return []
            if len(clade.imgids) == 1:
                cimg = list(clade.imgids)[0]
                if cimg in terminals:
                    #n = terminals[cimg]
                    return [(clname,1,1,
                                #n,n,n,
                                1.0)]
                return []

            core,pv,intersection = self.is_core( clade, terminals, er = er )

            if core:
                #ns = [terminals[ii] for ii in terminals_s if ii in clade.imgids]
                return [( clname,
                          len(intersection),len(clade.imgids),
                          #len(clade.imgids&terminals_s),len(clade.imgids),
                          #min(ns),max(ns),np.mean(ns),
                          pv)]
            rets = []
            for c in clade.clades:
                rets += _find_core_rec(c)
            return rets
        return  _find_core_rec( self.tree.root )


    def add_full_paths( self ):

        def _add_full_paths_( clade, path ):
            lpath = path + ([clade.name] if clade.name else [])
            clade.full_name = self.lev_sep.join( lpath )
            for c in clade.clades:
                _add_full_paths_( c, lpath )
        _add_full_paths_( self.tree.root, [] )

    def find_cores( self, cl_taxa_file, min_core_size = 1, error_rate = 0.95, subtree = None, skip_qm = True ):
        if subtree:
            self.subtree( 'name', subtree )
        self.ctc = {}
        imgids2terminals = {}
        for t in self.tree.get_terminals():
            #t.imgid = t.name # if "t__"in t.name else t.name
            t.imgid = t.name[3:] if "t__"in t.name else t.name # C2
            #t.imgid = int(t.name[3:] if "t__"in t.name else t.name) # C2
            t.nterminals = 1
            imgids2terminals[t.imgid] = t


        # can be made faster with recursion
        for n in self.tree.get_nonterminals():
            n.imgids = set( [nn.imgid for nn in n.get_terminals()]  )
            n.nterminals = len( n.imgids )

        self.add_full_paths() # unnecessary

        ret = {}
        for vec in (l.strip().split('\t') for l in open(cl_taxa_file)):
            sid = vec[0]
            #sid = int(vec[0]) # C2
            #tgts_l = [int(s) for s in vec[1:]]
            #tgts = dict([(s,tgts_l.count(s)) for s in set(tgts_l)])
            tgts = set([s for s in vec[1:]])
            #tgts = set([int(s) for s in vec[1:]]) C2

            if len(tgts) >= min_core_size:
                subtree_name = self.lev_sep.join(subtree.split(self.lev_sep)[:-1] ) if subtree else None
                ret[sid] = self._find_core( tgts, er = error_rate, root_name = subtree, skip_qm = skip_qm )
                #print sid #, ret[sid]
        return ret

    def markerness( self, coreness, uniqueness, cn_min, cn_max, cn_avg ):
        return coreness * uniqueness * (1.0 / float(cn_max-cn_min+1)) * 1.0 / cn_avg

    def find_markers( self, cu_file, hitmap_file, core_file ):
        self.ctc = {}
        imgids2terminals = {}
        ids2clades = {}
        for t in self.tree.get_terminals():
            t.imgid = t.name
            # t.imgid = int(t.name) # C2
            t.nterminals = 1
            imgids2terminals[t.imgid] = t
            ids2clades[t.name] = t

        # can be made faster with recursion (but it is not a bottleneck)
        for n in self.tree.get_nonterminals():
            n.imgids = set( [nn.imgid for nn in n.get_terminals()]  )
            n.nterminals = len( n.imgids )

        self.add_full_paths() # unnecessary

        #cus = dict([(int(l[0]),[int(ll) for ll in l[1:]]) for l in # C2
        cus = dict([(l[0],[ll for ll in l[1:]]) for l in
                        (line.strip().split('\t') for line in open(cu_file))])
        #cinfo = dict([(int(v[0]),[v[1]] + [int(vv) for vv in v[2:6]] + [float(vv) for vv in v[6:]]) # C2
        cinfo = dict([(v[0],[v[1]] + [vv for vv in v[2:6]] + [float(vv) for vv in v[6:]])
                        for v in (line.strip().split('\t') for line in open(core_file))])

        ret = {}
        for vec in (l.strip().split('\t') for l in open(hitmap_file)):
            sid = vec[0]
            #sid = int(vec[0]) # C2
            tgts_l = set([s for s in vec[1:]])
            # tgts_l = set([int(s) for s in vec[1:]]) # C2
            lca = self.lca( cus[sid], ids2clades )
            if lca.is_terminal():
                tin = set([lca.imgid])
                tout = tgts_l - tin
            else:
                tout = tgts_l - lca.imgids
                tin = lca.imgids & tgts_l
            ci = cinfo[sid]
            ltin = len(tin)
            ltout = len(tout)
            uniqueness = float(ltin)/float(ltin+ltout)
            coreness = float( ci[-1] )
            cn_min, cp_max, cn_avg = [float(f) for f in ci[-4:-1]]
            gtax = ci[0]
            cobs, ctot = ci[1], ci[2]
            # cobs, ctot = int(ci[1]), int(ci[2]) # C2
            markerness = self.markerness( coreness, uniqueness, cn_min, cp_max, cn_avg )

            res_lin = [ gtax, markerness, coreness, uniqueness, cobs, ctot, cn_min, cp_max, cn_avg,
                        ltin, ltout, "|".join([str(s) for s in tin]), "|".join([str(s) for s in tout]) ]
            ret[sid] = res_lin
        return ret


    def select_markers( self, marker_file, markerness_th = 0.0, max_markers = 200 ):
        cl2markers = colls.defaultdict( list )
        for line in (l.strip().split('\t') for l in open( marker_file )):
            gid = line[1]
            markerness = float(line[2])
            if markerness < markerness_th:
                continue
            cl2markers[gid].append( line )
        for k,v in cl2markers.items():
            cl2markers[k] = sorted(v,key=lambda x:float(x[2]),reverse=True)[:max_markers]
        return cl2markers.values()

    def get_c2t( self ):
        tc2t = {}

        def _get_c2t_( clade ):
            lterms = clade.get_terminals()
            tc2t[clade] = set([l.name for l in lterms])
            if clade.is_terminal():
                return
            for c in clade.clades:
                _get_c2t_( c )
        _get_c2t_( self.tree.root )
        return tc2t

    def ltcs( self, terminals, tc2t = None, terminals2clades = None, lca_precomputed = None ):
        set_terminals = set( terminals )
        lca = lca_precomputed if lca_precomputed else self.lca( terminals, terminals2clades )
        def _ltcs_rec_( clade, cur_max ):
            if clade.is_terminal() and clade.name in set_terminals:
                return clade,1
            terms = tc2t[clade] if tc2t else set([cc.name for cc in clade.get_terminals()])
            if len(terms) < cur_max:
                return None,0
            if terms <= set_terminals:
                return clade,len(terms)
            rets = []
            for c in clade.clades:
                r,tmax = _ltcs_rec_( c, cur_max )
                if tmax >= cur_max:
                    cur_max = tmax
                    if r:
                        rets.append((r,tmax))
            if rets:
                return sorted(rets,key=lambda x:x[1])[-1][0],cur_max
            else:
                return None,None
        return _ltcs_rec_( lca, cur_max = 0 )[0]

    def lca( self, terminals, terminals2clades = None ):
        clade_targets = []
        if terminals2clades:
            clade_targets = [terminals2clades[str(t)] for t in terminals]
        else:
            clade_targets = [t for t in self.tree.get_terminals() if t.name in terminals]
            """
            for t in terminals:
                ct = list(self.tree.find_clades( {"name": str(t)} ))
                if len( ct ) > 1:
                    sys.stderr.write( "Error: non-unique target specified." )
                    sys.exit(-1)
                clade_targets.append( ct[0] )
            """
        lca = self.tree.common_ancestor( clade_targets )
        return lca


    def lcca( self, t, t2c ):
        node_path = list(self.tree.get_path(t))
        if not node_path or len(node_path) < 2:
            return None,None,None
        tlevs = t2c[t].split(self.lev_sep)[2:-1]
        for p in node_path[-15:]:
            terms = list(p.get_terminals())
            descn = [t2c[l.name].split(self.lev_sep)[2:-1] for l in  terms if l.name!=t]
            if not descn or len(descn) < 2:
                continue

            l = tlevs[-1]
            descr_l = [d[-1] for d in descn]
            if len(set(descr_l)) == 1 and descr_l[0] != l and \
                l != "s__sp_" and not l.endswith("unclassified") and \
                descr_l[0] != "s__sp_" and not descr_l[0].endswith("unclassified"):
                return p,terms,self.lev_sep.join(tlevs)
        return None,None,None


    def tax_precision( self, c2t_f, strategy = 'lca' ):
        c2t = self.read_tax_clades( c2t_f )
        res = []
        for c,terms in c2t.items():
            lca = self.lca( terms )
            num = partial_branch_length(lca,terms)
            den = lca.total_branch_length()
            prec = num / den
            res.append([c,str(prec)])
        return res

    def tax_recall( self, c2t_f ):
        c2t = self.read_tax_clades( c2t_f )
        res = []
        for c,terms in c2t.items():
            lca = self.lca( terms )
            ltcs = self.ltcs( terms )
            lca_terms = set(lca.get_terminals())
            ltcs_terms = set(ltcs.get_terminals())
            out_terms = lca_terms - ltcs_terms
            outs = [c]
            if len(out_terms):
                diam = sum(sorted(ltcs.depths().values())[-2:])
                outs += [":".join([t.name,str( self.tree.distance(ltcs,t)/diam )])
                             for t in out_terms]
            res.append( outs )
        return res

    def tax_resolution( self, terminals ):
        pass


    def prune( self, strategy = 'lca', n = None, fn = None, name = None, newname = None ):
        prune = None
        if strategy == 'root_name':
            ct = list(self.tree.find_clades( {"name": name} ))
            if len( ct ) > 1:
                sys.stderr.write( "Error: non-unique target specified." )
                sys.exit(-1)
            prune = ct[0]
        elif strategy == 'lca':
            terms = self.read_targets( fn ) if isinstance(fn,str) else fn
            prune = self.lca( terms )
        elif strategy == 'ltcs':
            terms = self.read_targets( fn ) if isinstance(fn,str) else fn
            prune = self.ltcs( terms )
        elif strategy == 'n_anc':
            if n is None:
                n = 1
            ct = list(self.tree.find_clades( {"name": name} ))
            if len( ct ) > 1:
                sys.stderr.write( "Error: non-unique target specified.\n" )
                sys.exit(-1)
            node_path = list(self.tree.get_path(name))
            if not node_path or len(node_path) < n:
                sys.stderr.write( "Error: no anchestors or number of anchestors < n." )
                sys.exit(-1)
            toprune = node_path[-n]
            fat = node_path[-n-1]
            fat.clades = [cc for cc in fat.clades if cc != toprune]
            prune = None
        else:
            sys.stderr.write( strategy + " not supported yet." )
            sys.exit(-1)
        if prune:
            prune.clades = []
            if newname:
                prune.name = newname

    def subtree( self, strategy, n = None, fn = None ):
        newroot = None
        if strategy == 'name':
            ct = list(self.tree.find_clades( {"name": fn} ))
            if len( ct ) != 1:
                int_clades = self.tree.get_nonterminals()
                for cl in int_clades:
                    if n == cl.full_name:
                        ct = [cl]
                        break
                if not ct:
                    sys.stderr.write( "Error: target not found." )
                    sys.exit(-1)
            newroot = ct[0]
        elif strategy == 'lca':
            terms = self.read_targets( fn ) if isinstance(fn,str) else fn
            newroot = self.lca( terms )
        elif strategy == 'ltcs':
            terms = self.read_targets( fn ) if isinstance(fn,str) else fn
            newroot = self.ltcs( terms )
        if newroot:
            self.tree.root = newroot

    def rename( self, strategy, n = None, terms = None ):
        newroot = None
        if strategy == 'root_name':
            ct = list(self.tree.find_clades( {"name": n} ))
            if len( ct ) > 1:
                sys.stderr.write( "Error: non-unique target specified.\n" )
                sys.exit(-1)
            newroot = ct[0]
        elif strategy == 'lca':
            newroot = self.lca( terms )
        elif strategy == 'ltcs':
            newroot = self.ltcs( terms )
        if newroot:
            newroot.name = n

    def export( self, out_file ):
        self.tree = self.tree.as_phyloxml()
        Phylo.write( self.tree, out_file, "phyloxml")

    def read_tax_clades( self, tf ):
        with open( tf ) as inpf:
            return dict([(ll[0],ll[1:]) for ll in [l.strip().split('\t') for l in inpf]])

    def read_targets( self, tf ):
        if tf.count(":"):
            return tf.split(":")
        with open( tf ) as inpf:
            return [l.strip() for l in inpf]

    def reroot( self, strategy = 'lca', tf = None, n = None ):
        if strategy in [ 'lca', 'ltcs' ]:
            targets = self.read_targets( tf )

            lca = self.lca( targets ) if strategy == 'lca' else self.ltcs( targets )
            reroot_mid_fat_edge( self.tree, lca)

            #lca_f = get_parent( self.tree, lca )
            #
            #bl = lca.branch_length
            #new_clade = PClade(branch_length=bl*0.5, clades = [lca])
            #lca.branch_length = bl*0.5
            #if lca_f:
            #    lca_f.clades = [c for c in lca_f.clades if c != lca] + [new_clade]
            #    reroot( self.tree, new_clade )
            #else:
            #    self.tree.root = new_clade

        elif strategy == 'midpoint':
            pass
            #self.tree.reroot_at_midpoint(update_splits=True)
        elif strategy == 'longest_edge':
            nodes = list(self.tree.get_nonterminals()) + list(self.tree.get_terminals())
            longest = max( nodes, key=lambda x:x.branch_length )
            reroot_mid_fat_edge( self.tree, longest )
            #longest_edge = max( self.ntree.get_edge_set(),
            #                    key=lambda x:x.length)
            #self.tree.reroot_at_edge(longest_edge, update_splits=True)
        elif strategy == 'longest_internal_edge':
            nodes = list(self.tree.get_nonterminals())
            longest = max( nodes, key=lambda x:x.branch_length )
            if self.tree.root != longest:
                reroot_mid_fat_edge( self.tree, longest )
            #longest = get_lensorted_int_edges( self.tree )
            #self.tree.reroot_at_edge( list(longest)[-1], update_splits=True )
        elif strategy == 'longest_internal_edge_n':
            nodes = list(self.tree.get_nonterminals())
            longest = max( nodes, key=lambda x:x.branch_length
                    if len(x.get_terminals()) >= n else -1.0)
            reroot_mid_fat_edge( self.tree, longest )
            #longest = get_lensorted_int_edges( self.tree, n )
            #self.tree.reroot_at_edge( list(longest)[-1], update_splits=True )


    def reorder_tree( self ):
        self._ord_terms = []
        def reorder_tree_rec( clade ):
            if clade.is_terminal():
                self._ord_terms.append( clade )
                return clade,clade
            clade.clades.sort( key=lambda x:len(x.get_terminals()), reverse = True)
            for c in clade.clades:
                c.fc,c.lc = reorder_tree_rec( c )
            return clade.clades[0].fc,clade.clades[-1].lc
            #clade.fc, clade.lc = clade.clades[0], clade.clades[-1]

        reorder_tree_rec( self.tree.root )
        last = None # self._ord_terms[-1]
        for c in self._ord_terms:
            c.pc = last
            if last:
                last.nc = c
            last = c
        c.nc = None
        #self._ord_terms[-1].nc = None # self._ord_terms[0]


    def get_subtree_leaves( self, full_names = False ):

        subtrees = []
        def rec_subtree_leaves( clade ):
            if not len(clade.clades):
                return [clade.name]
            leaves = []
            for c in clade.clades:
                leaves += rec_subtree_leaves( c )
            leaves = [l for l in leaves if l]
            subtrees.append( (clade.name if clade.name else "",leaves)  )
            return leaves

        rec_subtree_leaves( self.tree.root )
        return subtrees

    def get_clade_names( self, full_names = False, leaves = True, internals = True ):
        clades = []
        if leaves:
            clades += self.tree.get_terminals()
        if internals:
            clades += self.tree.get_nonterminals()
        if full_names:
            def rec_name( clade, nam = "" ):
                ret = []
                if not nam and not clade.name:
                    lnam = ""
                elif not nam:
                    lnam = clade.name
                elif not clade.name:
                    lnam = nam
                else:
                    lnam = self.lev_sep.join( [nam, clade.name if clade.name else ""] )
                ret += [lnam] if lnam else []
                for c in clade.clades:
                    ret += rec_name(c,lnam)
                return ret
            names = set(rec_name(self.tree.root))
        else:
            names = set([c.name for c in clades])
        return sorted(names)

