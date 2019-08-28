import time
import argparse as ap
from utils import info

"""
Reads and parses the command line arguments of the script.

:returns: the parsed arguments
"""
def read_params():
    p = ap.ArgumentParser(description="")
    p.add_argument('-i', '--input', type=str, default=None,
                   help="The input folder with the read files")
    
    return p.parse_args()


def parse_psortb(input):
    c=0
    with open(input, 'rt') as ifn:
        for line in ifn:
            line = line.strip().split("\t")
            # info(line[1], init_new_line=True)
            if line[1] == "Periplasmic":
                c+=1
                print("cp ../../04_coli/s__Escherichia_coli/"+line[0].split()[0].replace("s--Escherichia-coli_UniRef90-","").replace("_0", "")+".faa .")
        info(c, init_new_line=True)

"""
Main call

:param input: the file to parse
"""
if __name__ == "__main__":
    t0 = time.time()
    args = read_params()
    # info(format(time.ctime(int(time.time())))+" Start execution")

    parse_psortb(args.input)

    exec_time = time.time() - t0
    info(format(time.ctime(int(time.time())))+" Finish execution ("+str(round(exec_time, 2))+
        " seconds)\n", init_new_line=True)