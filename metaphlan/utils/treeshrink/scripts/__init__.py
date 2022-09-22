#!/usr/bin/env python

#############################################################################
##  this file is part of TreeShrink.
##  see "license.txt" for terms and conditions of usage.
#############################################################################

PROGRAM_NAME = "TREESHRINK"
PROGRAM_AUTHOR = ["Uyen Mai","Siavash Mirarab"]
PROGRAM_LICENSE = "GNU General Public License, version 3"
PROGRAM_VERSION = "1.3.9"
PROGRAM_YEAR = "2017"
PROGRAM_DESCRIPTION = "Fast and accurate detection of outlier long branches in collections of phylogenetic trees"
PROGRAM_WEBSITE = "https://uym2.github.io/TreeShrink/"
PROGRAM_INSTITUTE = "Department of Computer Science and Engineering, University of California at San Diego"


from tempfile import mkdtemp,mktemp
from os.path import basename, dirname, splitext,realpath,join,normpath,isdir,isfile,exists
from os import mkdir,rmdir
        
tempdir = None
def get_tmp_dir():
    global tempdir
    return tempdir
def set_tmp_dir(t):
    global tempdir
    if t:
        tempdir = t
        try:
            mkdir(tempdir)
        except FileExistsError:
            pass
    else:
        tempdir = mkdtemp()
    return tempdir
def get_tmp_file(name=None,prefix=None):
    global tempdir
    if name is None:
        if prefix:
            return mktemp(dir=tempdir,prefix=prefix)
        else:
            return mktemp(dir=tempdir)
    else:
        return normpath(join(tempdir,name))
        
    

import sys
sys.setrecursionlimit(5000)
