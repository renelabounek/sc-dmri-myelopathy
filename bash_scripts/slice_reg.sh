#!/bin/bash

#============================================================================================
# Script for runnning translation registration for correction motion corrupted images.
# This script is called from matlab script T2TRA_shift_correction_interp.m
# JV, 2018
#============================================================================================
#
# Copyright 2016-2020 Rene Labounek (1,2,3,4), Jan Valosek (1,2) and Petr Hlustik (1,2)
#
# 1 - University Hospital Olomouc, Olomouc, Czech Republic
# 2 - Palacky University Olomouc, Olomouc, Czech Republic
# 3 - University Hospital Brno, Brno, Czech Republic 
# 4 - University of Minnesota, Minneapolis, US
#
# This file is part of sc-dmri-myelopathy available at: https://github.com/renelabounek/sc-dmri-myelopathy
#
# Please, cite sc-dmri-myelopathy as:
# Labounek R, Valosek J, Horak T, Svatkova A, Bednarik P, Vojtisek L, Horakova M, Nestrasil I,
# Lenglet C, Cohen-Adad J, Bednarik J and Hlustik P. HARDI-ZOOMit protocol improves specificity
# to microstructural changes in presymptomatic myelopathy. Scientific Reports [Revised; Under review]
#
# sc-dmri-myelopathy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# sc-dmri-myelopathy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with sc-dmri-myelopathy.  If not, see <https://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------------------------------------

SCTDIR=/usr/local/lib/sct/bin
PATH=$PATH:$SCTDIR


fixed=$1
moving=$2

echo "Starting slice-by-slice registration by isct_antsSliceRegularizedRegistration function..."
isct_antsSliceRegularizedRegistration -t Translation[0.5] -m CC[$fixed,$moving,1,4,Regular,0.2] -p 5 -i 10 -f 1 -s 5 -v 0 -o [slice_reg,odd_to_even.nii] > isct_antsSliceRegularizedRegistration.log 2>&1
#echo $fixed $moving
echo "Slice-by-slice registration by isct_antsSliceRegularizedRegistration function is completed."
