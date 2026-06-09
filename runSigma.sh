#!/bin/bash

# source /full_path_to/bin/thisroot.sh # if needed
root -b -l -x <<EOF
.L Sigma.C
Sigma a;
a.Loop();

EOF
