#!/bin/bash

#============================================================================================
# Script for runnning translation registration for correction motion corrupted images.
# This script is called from matlab script T2TRA_shift_correction_interp.m
# JV, 2018
#============================================================================================

SCTDIR=/usr/local/lib/sct/bin
PATH=$PATH:$SCTDIR


fixed=$1
moving=$2

echo "Starting slice-by-slice registration by isct_antsSliceRegularizedRegistration function..."
isct_antsSliceRegularizedRegistration -t Translation[0.5] -m CC[$fixed,$moving,1,4,Regular,0.2] -p 5 -i 10 -f 1 -s 5 -v 0 -o [slice_reg,odd_to_even.nii] > isct_antsSliceRegularizedRegistration.log 2>&1
#echo $fixed $moving
echo "Slice-by-slice registration by isct_antsSliceRegularizedRegistration function is completed."
