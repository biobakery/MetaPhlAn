#!/usr/bin/env python

import sys
import os
import bz2
import gzip

p2= float(sys.version_info[0]) < 3.0

def ignore_spaces( l ):
    if l[0] == '@' or l[0] == '>':
        return l.replace(' ','_')
    return l

def fopen( fn ):
    fileName, fileExtension = os.path.splitext( fn )
    if fileExtension == '.bz2':
        return (bz2.open(fn, "rt") if not p2 else bz2.BZ2File(fn,"r"))
    if fileExtension == '.gz':
        return gzip.open(fn,"rt") 
    return open(fn)


if __name__ == '__main__':

    if len(sys.argv) > 1:
        if ((sys.argv[1] == '-h') or (sys.argv[1] == '--help') or
            (sys.argv[1] == '-v') or (sys.argv[1] == '--version')):
            sys.exit()
    
    if len(sys.argv) < 2:
        for l in sys.stdin:
            sys.stdout.write(ignore_spaces(l))
    
    if len(sys.argv) > 1:
        for a in sys.argv[1:]:
            for f in a.split(','):
                with fopen(f) as inf:
                    for l in inf:
                        sys.stdout.write(ignore_spaces(l))

