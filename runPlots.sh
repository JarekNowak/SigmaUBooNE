#!/bin/bash

# source /full_path_to/bin/thisroot.sh # if needed

# Generates the stacked + shape PNGs from SigmaHists.root (produced by runSigma.sh).
# Compiled with ACLiC (the trailing '+') for speed and compile-time error checking.
root -b -l -q PlotSigma.C+
