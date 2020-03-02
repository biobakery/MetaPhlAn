#!/bin/bash


# VTML240 60 (MUSCLE)
perl vt_scores.pl VTML 240 vtml240 60
python3 serialize_matrix.py vtml240

# VTML 200 3
perl vt_scores.pl VTML 200 vtml200
python3 serialize_matrix.py vtml200

# MIQS
python3 serialize_matrix.py miqs

# PFASUM60
python3 serialize_matrix.py pfasum60
