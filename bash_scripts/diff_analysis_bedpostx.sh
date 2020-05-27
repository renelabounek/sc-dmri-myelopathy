#!/bin/bash

#-----------------------------------------------------------------------------------------------------------
# Script for estimation of Ball-and-Sticks model on the spinal cord diffusion data (ZOOMit and/or RESOLVE)
# USAGE:
#			diff_analysis_bedpostx.sh <DATA_FOLDER> <SUB_ID> <SEQ_order>
#						- DATA_FOLDER - folder containing individual subjects data
#						- SUB_ID - name of the subject folder with raw (nifti) anatomical and diffusion data
#						- SEQORDER - order of diffusion protocols (ZOOMit_interp, ZOOMit_interp_moco, ZOOMit_nointerp, ZOOMit_nointerp_moco, RESOLVE)
# EXAMPLE:
#			diff_analysis_bedpostx.sh /home/user/data 1234B 11001
#
# Analysis contains:
#		1) Estimation of Ball-and-Sticks model using FSL's bedpostX function
#		2) Warping of estimated maps/metrics into anatomical T2 transversal (T2TRA) space.
#
# Jan Valosek, Rene Labounek. fMRI laboratory Olomouc, Czech Republic, 2017-2020
#
#-----------------------------------------------------------------------------------------------------------
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

# Some global variable
VERSION=01-03-2020

SCT_VER_GIT="7af1436782a5ea87758ac893a5f658d58eb3c60a"		# SCT v3.2.3 - https://github.com/neuropoly/spinalcordtoolbox/commit/7af1436782a5ea87758ac893a5f658d58eb3c60a

SCRIPT_NAME=${0##*/}													# get script name
SCRIPT_DIR="$(cd "$(dirname $0)";pwd -P)"			# get directory where is this script located
LINE="========================================================================================================================================="

BINDIR=${SCRIPT_DIR}												# our bash scripts

pidlist=""		# list of process IDs

# Trap kill from check_input function if SIGUSR1 is set
trap "echo Exitting>&2;exit" SIGUSR1

# Initialization and print help
init()
{

	# Fetch some basic functions
	source ${BINDIR}/bash_basic_functions.sh

	if [[ $# -eq 0 ]] || [[ $1 == "-h" ]];then
		echo -e "$LINE\nHelp for $SCRIPT_NAME script performing estimation of Ball-and-sticks model (using bedpostX FSL function)"
		echo -e "and warping estimated metrics/maps from diffusion to anatomical space, version $VERSION"
		echo -e "REQUIREMENTS: Installed bash interpreter, FSL and Spinal Cord Toolbox libraries and HTCondor.\n$LINE"
		echo -e "USAGE:\n\t$SCRIPT_NAME <DATA_FOLDER> <SUB_ID> <SEQ_order>\nEXAMPLE:\n\t$SCRIPT_NAME /home/user/data 1234B 11001"
		echo -e "or\n\t$SCRIPT_NAME/md2/NA-CSD subjects.txt\n$LINE"
		echo -e "It is neccesary to run this script on condor master machine\n$LINE"
		echo -e "Valosek, Labounek 2017-2020\tfMRI laboratory, Olomouc, CZ\n$LINE"
	elif [[ $# -eq 1 ]] && [[ $1 != "-h" ]] || [[ $# -gt 3 ]];then
		show "Invalid form of input argument!\nExiting!" et
	else

		Diff_preproc_ZOOMit_interp=${3:0:1}
		Diff_preproc_ZOOMit_interp_moco=${3:1:1}
		Diff_preproc_ZOOMit_notinterp=${3:2:1}
		Diff_preproc_ZOOMit_notinterp_moco=${3:3:1}
		RESOLVE_RUN=${3:4:1}

		ACQ=""
		if [[ $Diff_preproc_ZOOMit_interp == 1 ]];then
        		ACQ="$ACQ Diff_preproc_ZOOMit_interp"
		fi
		if [[ $Diff_preproc_ZOOMit_interp_moco == 1 ]];then
        		ACQ="$ACQ Diff_preproc_ZOOMit_interp_moco"
		fi
		if [[ $Diff_preproc_ZOOMit_notinterp == 1 ]];then
        		ACQ="$ACQ Diff_preproc_ZOOMit_notinterp"
		fi
		if [[ $Diff_preproc_ZOOMit_notinterp_moco == 1 ]];then
		        ACQ="$ACQ Diff_preproc_ZOOMit_notinterp_moco"
		fi


		DATA=$1
		if [[ ${DATA: -1} == "/" ]];then
			DATA=${DATA%/}			# delete slash at the end of DATA varaible if is containted
		fi

		check_input d ${DATA}

		isFile=$(file ${DATA}/$2 | cut -d\  -f2)		# Check if the second argument is txt file contating list of subjects
		if [[ $isFile = "ASCII" ]];then
			subjects=$(cat ${DATA}/$2)
		else
			subjects=$2
		fi

		for SUB in $subjects;do
			if [[ ${SUB: -1} == "/" ]];then
				SUB=${SUB%/}			# delete slash at the end of SUB varaible if is containted
			fi

			check_input d ${DATA}/${SUB}

			main ${DATA} ${SUB}
		done

	fi
}
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Main function
main()
{

		DATA=$1
		SUB=$2

		export LOGPATH=${DATA}/${SUB}/`echo $SCRIPT_NAME | sed 's/.\{3\}$//'`.log
		exec > >(tee -a $LOGPATH) 2>&1

		# Get SCT version
		SCT_VER=$(sct_check_dependencies | (grep "Spinal Cord Toolbox") | awk -F '/' '{print $NF}' | sed 's/)//')		# run sct_check_dependencies script, get line containing 'Spinal Cord Toolbox', delete everything before last slash and delete last bracket

		if [[ $SCT_VER != $SCT_VER_GIT ]]; then
			show "SCT v3.2.3 was not found. Current version of SCT is different - $SCT_VER\nDownload SCT v3.2.3 by:\ngit clone --branch v3.2.3 https://github.com/neuropoly/spinalcordtoolbox.git" e
		fi
		
		show "$LINE" yt
		show "\tStarting: $SCRIPT_NAME script, version: $VERSION" y
		show "\tData directory: ${DATA}" y
		show "\tSubject: ${SUB}" y
		show "\tMachine name: $HOSTNAME" y
		show "\tNum of CPU: $(nproc)" y
		show "\tSCT version: $SCT_VER" y
		show "$LINE" yt

		check_input b matlab fslmerge sct_dmri_moco

	#----------------------------------------------
	# bedpostx
		#ZOOMit
		NAME="zoomit"

		for FOLDER in $ACQ;do
			if [[ ! -d ${DATA}/${SUB}/Results/${FOLDER}/dMRI.bedpostX ]]; then
				show "Starting bedpostX in ${FOLDER}" y
				bedpostx_function ${DATA} ${FOLDER} $NAME > ${DATA}/${SUB}/Results/${FOLDER}/bedpostx.log &		# Call function for bedpostx_function
				pidlist="$pidlist $!"    # add process ID of most recent executed background pipeline
			else
				show "Bedpostx in ${FOLDER} folder has been done before!" g
			fi
		done

		if [[ $RESOLVE_RUN == 1 ]];then
			#RESOLVE
			NAME="resolve"
			FOLDER="Diff_preproc_RESOLVE"
			if [[ ! -d ${DATA}/${SUB}/Results/${FOLDER}/dMRI.bedpostX ]]; then
			show "Starting bedpostX in ${FOLDER}" y
				bedpostx_function ${DATA} ${FOLDER} $NAME > ${DATA}/${SUB}/Results/${FOLDER}/bedpostx.log &		# Call function for bedpostx_function
				pidlist="$pidlist $!"    # add process ID of most recent executed background pipeline
			else
				show "BedpostX in ${FOLDER} folder has been done before!" g
			fi
		fi

		if [[ "$pidlist" != "" ]]; then
			wait $pidlist
			pidlist=""
		fi

	#----------------------------------------------
	# masking
		#ZOOMit
		NAME="zoomit"

		for FOLDER in $ACQ;do
			if [[ ! -f ${DATA}/${SUB}/Results/${FOLDER}/dMRI.bedpostX/mean_f1samples_cord.nii* ]]; then
				show "Starting masking of bedpostX results in ${FOLDER}" y
				masking ${DATA} ${FOLDER} $NAME > ${DATA}/${SUB}/Results/${FOLDER}/masking.log &		# Call function for masking
				pidlist="$pidlist $!"    # add process ID of most recent executed background pipeline
			else
				show "Masking in ${FOLDER} folder has been done before!" g
			fi
		done

		if [[ $RESOLVE_RUN == 1 ]];then
			#RESOLVE
			NAME="resolve"
			FOLDER="Diff_preproc_RESOLVE"
			if [[ ! -f ${DATA}/${SUB}/Results/${FOLDER}/dMRI.bedpostX/mean_f1samples_cord.nii* ]]; then
				show "Starting masking of bedpostX results in ${FOLDER}" y
				masking ${DATA} ${FOLDER} $NAME > ${DATA}/${SUB}/Results/${FOLDER}/masking.log &		# Call function for masking
				pidlist="$pidlist $!"    # add process ID of most recent executed background pipeline
			else
				show "Masking in ${FOLDER} folder has been done before!" g
			fi
		fi

		if [[ "$pidlist" != "" ]]; then
			wait $pidlist
			pidlist=""
		fi

	#----------------------------------------------
		show "$LINE" yt
		show "\t$SCRIPT_NAME for ${DATA}/${SUB} finished." y
		show "$LINE" yt

		# END of main function
}

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Function which is called for estimation Ball-and-Sticks model using bedpostX FSL function
bedpostx_function()
{
	DATA=$1
	FOLDER=$2
	NAME=$3

	cd ${DATA}/${SUB}/Results/${FOLDER}/

	if [[ ! -d ${DATA}/${SUB}/Results/${FOLDER}/dMRI ]];then
		show "${DATA}/${SUB}/Results/${FOLDER}/dMRI does not exist! Run diff_analysis.sh first." e
	fi

	show "Starting bedpostX in ${DATA}/${SUB}/Results/${FOLDER} folder!" y

	bedpostx_condor dMRI -n 2 -model 1 > ${DATA}/${SUB}/Results/${FOLDER}/bedpostx_condor_dagman.log
	sleep 30
	#JOBID=`sed '6q;d' ${DATA}/${SUB}/Results/${FOLDER}/bedpostx_condor_dagman.log | sed 's/.* //' | sed 's/.$//'`
	while [[ ! -f ${DATA}/${SUB}/Results/${FOLDER}/dMRI.bedpostX/condor_logs/bedpostx_postproc_condor.out ]];do
		sleep 60
	done
	sleep 150

	show "BedpostX in ${FOLDER} folder is done!" y
}

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Masking of results from Ball-and-Sticks model and warping these results into T2TRA space
masking()
{
	DATA=$1
	FOLDER=$2
	NAME=$3

	if [[ ! -f ${DATA}/${SUB}/Results/${FOLDER}/dMRI.bedpostx/dyads1.nii.gz ]];then
		show "Starting postprocessing of bedpostx data in ${DATA}/${SUB}/Results/${FOLDER} folder!" y
		${DATA}/${SUB}/Results/${FOLDER}/dMRI.bedpostX/condor_logs/bedpostx_postproc_condor.sh
	fi

	show "Starting masking of diffusion data in ${DATA}/${SUB}/Results/${FOLDER} folder!" y

	cd ${DATA}/${SUB}/Results/${FOLDER}

	fslmaths dMRI.bedpostX/mean_f1samples.nii.gz -mas dwi_mean_seg.nii.gz dMRI.bedpostX/mean_f1samples_cord.nii.gz
	#gunzip dMRI.bedpostX/mean_f1samples_cord.nii.gz
	fslmaths dMRI.bedpostX/mean_f2samples.nii.gz -mas dwi_mean_seg.nii.gz dMRI.bedpostX/mean_f2samples_cord.nii.gz
	#gunzip dMRI.bedpostX/mean_f2samples_cord.nii.gz
	fslmaths dMRI.bedpostX/dyads2_thr0.05.nii.gz -mas dwi_mean_seg.nii.gz dMRI.bedpostX/dyads2_thr0.05_cord.nii.gz
	#gunzip dMRI.bedpostX/dyads2_thr0.05_cord.nii.gz

	show "masking of diffusion data in ${FOLDER} folder is done!" y

	# Warping DIFF results into T2TRA
	show "Starting warping of diffusion data into T2TRA in ${DATA}/${SUB}/Results/${FOLDER} folder!" y

	check_input dc ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra

	for FILE in dMRI.bedpostX/mean_f1samples.nii.gz dMRI.bedpostX/mean_f1samples_cord.nii.gz dMRI.bedpostX/mean_f2samples.nii.gz dMRI.bedpostX/mean_f2samples_cord.nii.gz dMRI.bedpostX/dyads2_thr0.05.nii.gz dMRI.bedpostX/dyads2_thr0.05_cord.nii.gz dMRI.bedpostX/dyads1.nii.gz dMRI.bedpostX/mean_dsamples.nii.gz;do
		OUTPUT=`echo $FILE | sed 's:.*/::g'`
		sct_apply_transfo -i ${DATA}/${SUB}/Results/${FOLDER}/${FILE} -d ${DATA}/${SUB}/Results/Anat_Preproc/T2TRA_thr_bias_corr.nii.gz -w ${DATA}/${SUB}/Results/${FOLDER}/warp_DIFFtoT2TRA.nii.gz -o ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra/$OUTPUT
	done

	fslmaths ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra/mean_f1samples_cord.nii.gz -thr 0.1 ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra/mean_f1samples_cord.nii.gz
	fslmaths ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra/mean_f2samples_cord.nii.gz -thr 0.05 ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra/mean_f2samples_cord.nii.gz
	fslmaths ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra/dyads2_thr0.05_cord.nii.gz -mas ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra/mean_f2samples_cord.nii.gz ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra/dyads2_thr0.05_cord.nii.gz

	show "Warping of diffusion data into T2TRA in ${DATA}/${SUB}/Results/${FOLDER} folder is done!" y

}
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

init $@
