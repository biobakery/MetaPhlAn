#! /bin/bash


# remove anaconda from PATH
export PATH=$HOME/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
export PATH=/shares/CIBIO-Storage/CM/mir/tools/bin:$PATH

python metaphlan2.py $*
