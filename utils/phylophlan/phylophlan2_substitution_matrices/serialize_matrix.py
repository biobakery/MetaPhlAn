#!/usrbin/env python3


import pickle
import os
import sys


mat = {}
header = []

if not os.path.isfile(sys.argv[1]+'.pkl'):
    with open(sys.argv[1]) as f:
        for row in f:
            if not row.startswith('#'):
                row_clean = row.strip().split()

                if not header:
                    header = row_clean
                else:
                    aa = row_clean[0]

                    for i, s in enumerate(row_clean[1:]):
                        mat[(header[i], aa)] = float(s)

    with open(sys.argv[1]+'.pkl', 'wb') as f:
        pickle.dump(mat, f, protocol=pickle.HIGHEST_PROTOCOL)
