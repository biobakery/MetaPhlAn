# I/O lib for molecular sequences

from os.path import isfile
from os import remove
from random import random
from copy import copy
from . import get_tmp_file

try:
    import cPickle as pickle
except:
    import pickle


def get_taxon_list(filename):
    taxon_list = []
    for line in open(filename,'r'):
        if line[0] == '>':
            taxon_list = taxon_list + [line[1:].rstrip()]
    return sorted(taxon_list)


def hash_taxon_seq(filename):
    taxon_dict = {}
    f = open(filename,'r')
    for line in f:
        if line[0] == '>':
            taxon_dict[line[1:-1]] = gap_rm(f.next().rstrip())
    return taxon_dict

def gap_rm(str0,gap='-'):
    str1 = ''
    for c in str0:
        if c != gap:
            str1 =  str1 + c
    return str1


indexfiles = {}
def get_index_file_name(file_in):
    try:
        return indexfiles [file_in]
    except KeyError:
        tmp = get_tmp_file()
        indexfiles[file_in] = tmp
        return tmp

def index_fasta(file_in,file_out=None,store_index_file=True):
    # only work for fasta format
    f = open(file_in,'r')
    seq_pointers = {}
    count = {}
    fp = 0
    while 1:
        line = f.readline()
        if not line:
            break
        if line[0] == '>':
            seqName = line.rstrip()[1:]
            if seqName in count:
                c = count[seqName]
                p = float(c)/(c+1)
                r = random()
                #print(r)
                count[seqName] += 1   
                if r <= p:
                    continue    
            else:
                count[seqName] = 1
                    
            seq_pointers[seqName] = fp
        fp = f.tell()

    if not store_index_file:
        return seq_pointers
    
    if not file_out:
        file_out = get_index_file_name(file_in)
    with open(file_out,'wb') as fout:
        pickle.dump(seq_pointers,fout)
    f.close()

    return seq_pointers

def load_index(file_in,store_index_file=True,renew_index_file=False):
    file_idx = get_index_file_name(file_in)

    if renew_index_file or not isfile(file_idx):
        if renew_index_file:
            remove(file_idx)
        seq_pointers = index_fasta(file_in,store_index_file=store_index_file)
    else:
        with open(file_idx,'rb') as f:
            seq_pointers = pickle.load(f)

    return seq_pointers

def sample_from_list(file_in,taxa_list,file_out,store_index_file=True,renew_index_file=False):
    seq_pointers = load_index(file_in,store_index_file=store_index_file,renew_index_file=renew_index_file)
    with open(file_in,'r') as fin:
        with open(file_out,'w') as fout:     
            for taxon in taxa_list:
                try:
                    fin.seek(seq_pointers[taxon])
                    fout.write(fin.readline())
                    L = fin.readline()
                    while L[0] != '>':
                        fout.write(L.rstrip())
                        L = fin.readline()
                        if not L:
                            break
                    fout.write('\n')        
                except KeyError as e:
                    print ('taxon inconsistent in alignment and tree files: %s' %file_in )
                    raise e

def filter_out_by_list(file_in,removing_list,file_out,store_index_file=True,renew_index_file=False):
    seq_pointers = load_index(file_in,store_index_file=store_index_file,renew_index_file=renew_index_file)
    taxa_list = list(set(seq_pointers.keys()) - set(removing_list))
    
    with open(file_in,'r') as fin:
        with open(file_out,'w') as fout:     
            for taxon in taxa_list:
                try:
                    fin.seek(seq_pointers[taxon])
                    fout.write(fin.readline())
                    fout.write(fin.readline())
                except KeyError as e:
                    print ('taxon inconsistent in alignment and tree files: %s' %file_in)
                    raise e


def count_gaps(seq_aln):
    N = len(seq_aln[0])
    gap_count = [0]*N
    for seq in seq_aln:
        for i in range(N):
            gap_count[i] += (seq[i] == '-')
    return gap_count

def read_fasta(fas_file):
    taxon_names = []
    seq_aln = []
    is_first_seq = True
    with open(fas_file,'r') as f:
        for line in f:
            if line[0] == '>':
                taxon_names.append(line[1:].rstrip())
                if not is_first_seq:
                    seq_aln.append(new_seq)
                is_first_seq = False
                new_seq = ''
            else:
                new_seq += line.upper().rstrip()
    seq_aln.append(new_seq)
    return taxon_names, seq_aln    

def sort_aln(taxon_names,seq_aln):
    first_nongap_pos = [ len(seq_aln[0]) for i in range(len(seq_aln)) ]
    
    for i,seq in enumerate(seq_aln):
        for j,r in enumerate(seq):
            if r != '-':
                first_nongap_pos[i] = j
                break
    sorted_idx = sorted(range(len(first_nongap_pos)),key=lambda x:first_nongap_pos[x])
    sorted_names = []
    sorted_aln = []
    for i in sorted_idx:
        sorted_aln.append(seq_aln[i])
        sorted_names.append(taxon_names[i])

    return sorted_names, sorted_aln
    

def write_fasta(output_file,taxon_names,seq_aln):
    with open(output_file,'w') as f:
        T = len(taxon_names)
        for i in range(T):
            f.write(">"+taxon_names[i]+"\n")
            f.write(seq_aln[i]+"\n")

def is_aligned(fas_file):
    taxa,seqs = read_fasta(fas_file)
    l = len(seqs[0])
    for seq in seqs:
        if len(seq) != l:
            return False
    return True

def gap_propagate(cons_seq,targ_seq):
# propagate gaps from cons_seq to targ_seq; this is used after the cons_seq was aligned with another sequence and we want to
# propage the alignment to targ_seq. This is a similar idea with transitivity used in PASTA; the ultimate goal is to merge 2 alignments,
# but in this case we have a consensus sequence for each alignments and we also have a good way (properly using seconday structure) to align them.
# NOTE: gap_propagate is NOT symmetric: only propagate from consensus to target, not in reverse; careful consider which sequence is the cons_seq and which is targ_seq!

    out_seq = ''
    i = 0
    for c in cons_seq:
        if c != '-':
            out_seq += targ_seq[i]
            i += 1
        else:
            out_seq += '-'

    return out_seq

def impose_struct(pri_seq,str_seq):
    out_pri = ''
    out_str = ''

    for i,c in enumerate(pri_seq):
        if c != '-':
            out_pri += c.upper()
            c1 = str_seq[i]
            if c1 == '(' or c1 == '<' or c1 == '{':
                out_str += '('
            elif c1 == ')' or c1 == '>' or c1 == '}':
                out_str += ')'
            else:
                out_str += '.'
    
    return out_pri, out_str

def p_distance(seq1,seq2):
    d = 0
    count = len(seq1)
    for i,x in enumerate(seq1):
        y = seq2[i]
        if x == '-' and y == '-':
            count -= 1
        elif x != y:
            d = d+1
    return float(d)/count

def replace(from_letter, to_letter, aln):
    # replace all of the 'from_letter' to the 'to_letter' in aln
    locations = []
    new_aln = []
    for i,s in enumerate(aln):
        j_nongap = 0
        new_s = ''
        for x in s:
            if x == from_letter:
                new_s += to_letter 
                locations.append((i,j_nongap))
            else:
                new_s += x    
            j_nongap += (x != '-')
        new_aln.append(new_s)
    return new_aln,locations

def replace_back(letter, aln, locations):
    # replace the positions listed in 'locations' with the specified letter
    for i,j_nongap in locations:
        new_s = ''
        for x in aln[i]:
            if j_nongap == 0:
                new_s += letter
            else:
                new_s += x    
            j_nongap -= (x != '-')
        aln[i] = new_s        
            
def merge_rep_locations(rep_locations1,len1,rep_locations2):
    rep_locations = copy(rep_locations1)

    for i,j in rep_locations2:
        rep_locations.append((i+len1,j))

    return rep_locations
