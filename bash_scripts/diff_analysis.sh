#!/bin/bash

#-----------------------------------------------------------------------------------------------------------
# This script is called automatically from run_analysis.sh script or can be run manually.
# USAGE:
#			diff_analysis.sh <DATA_FOLDER> <SUB_ID> <SEQ_order>
#						- DATA_FOLDER - folder containing individual subjects data
#						- SUB_ID - name of the subject folder with raw (nifti) anatomical and diffusion data
#						- SEQORDER - order of diffusion protocols (ZOOMit_interp, ZOOMit_interp_moco, ZOOMit_nointerp, ZOOMit_nointerp_moco, RESOLVE)
# EXAMPLE:
#			diff_analysis.sh /home/user/data 1234B 11001
#
# Script for performing analysis of spinal cord diffusion data (ZOOMit and/or RESOLVE) only in diff space.
# Analysis contains:
#		1) Diffusion preprocessing (merge AP and PA b0 images) with or without motion correction by sct_dmri_moco
#		2) Correction of artifacts (geometrical distorsions and eddy current artifacts using FSL topup and eddy  functions (with whole FOV or with manually segmented mask of SC from topup_mean image))
#		3) Estimation of DTI model (using FSL function dtifit)
#		4) Registration between T2TRA and DIFF spaces - 2-step registration (see description of this function)
#		5) vertebrae labeling in DIFF space
#		6) masking results from dtifit by mask of SC
#
#	!!! It is necessary to run run_analysis.sh before this script to gain T2TRA_thr_bias_corr_seg.nii.gz !!!
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
# the Free Software Foundation, either version 3 of the License, or any later version.
#
# sc-dmri-myelopathy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with sc-dmri-myelopathy.  If not, see <https://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------------------------------------

# Some global variables
VERSION=01-03-2020

SCT_VER_GIT="7af1436782a5ea87758ac893a5f658d58eb3c60a"		# SCT v3.2.3 - https://github.com/neuropoly/spinalcordtoolbox/commit/7af1436782a5ea87758ac893a5f658d58eb3c60a

SCRIPT_NAME=${0##*/}													# get script name
SCRIPT_DIR="$(cd "$(dirname $0)";pwd -P)"			# get directory where is this script located
LINE="===================================================================================================================="

DEEPSEG=0			# Type of SC segmentation: 0 - sct_propseg, 1 - sct_deepseg_sc

BINDIR=${SCRIPT_DIR}												# our bash scripts

pidlist=""		# list of process IDs

# Trap kill from check_input function if SIGUSR1 is set
trap "echo Exiting>&2;exit" SIGUSR1

# Initialization and print help
init()
{
		# Fetch some basic functions
		source ${BINDIR}/bash_basic_functions.sh

		if [ $# -eq 0 ] || [ $1 == "-h" ];then
			echo -e "$LINE\nHelp for script performing analysis of the spinal cord diffusion data (ZOOMit and/or RESOLVE), version $VERSION."
			echo -e "Analysis contains:\n\tdiffusion preprocessing (merge AP and PA b0 images) with or withnout motion correction by sct_dmri_moco\n\ttopup and eddy (with whole FOV or with manually segmented mask of SC from topup_mean image)\n\tDTI estimation (using dtifit)\n\tregistration between T2TRA and DIFF spaces\n\tvertebrae labeling in DIFF space\n\tmasking results from dtifit by mask of SC"
			echo -e "REQUIREMENTS: Installed bash interpreter, Matlab, FSL and Spinal Cord Toolbox libraries.\n$LINE"
			echo -e "USAGE:\n\n$SCRIPT_NAME <DATA_FOLDER> <SUB_ID> <SEQ_order>\n\nEXAMPLE:\n\n$SCRIPT_NAME /home/user/data 1234B 11001"
			echo -e "or\n$SCRIPT_NAME /md2/NA-CSD subjects.txt\n$LINE"
			echo -e "Valosek, Labounek 2018\tfMRI laboratory, Olomouc, CZ\n$LINE"
		elif [ $# -eq 1 ] && [ $1 != "-h" ] || [ $# -gt 3 ];then
			check_input e "Invalid form of input argument(s)!"
		else

		Diff_preproc_ZOOMit_interp=${3:0:1}
		Diff_preproc_ZOOMit_interp_moco=${3:1:1}
		Diff_preproc_ZOOMit_notinterp=${3:2:1}
		Diff_preproc_ZOOMit_notinterp_moco=${3:3:1}
		RESOLVE_RUN=${3:4:1}

		ACQ=""
		if [ $Diff_preproc_ZOOMit_interp == 1 ];then
	        	ACQ="$ACQ Diff_preproc_ZOOMit_interp"
		fi
		if [ $Diff_preproc_ZOOMit_interp_moco == 1 ];then
	        	ACQ="$ACQ Diff_preproc_ZOOMit_interp_moco"
		fi
		if [ $Diff_preproc_ZOOMit_notinterp == 1 ];then
	        	ACQ="$ACQ Diff_preproc_ZOOMit_notinterp"
		fi
		if [ $Diff_preproc_ZOOMit_notinterp_moco == 1 ];then
	        	ACQ="$ACQ Diff_preproc_ZOOMit_notinterp_moco"
		fi

		DATA=$1
		if [ ${DATA: -1} == "/" ];then
			DATA=${DATA%/}			# delete slash at the end of DATA varaible if is containted
		fi

		check_input d $DATA

		isFile=$(file $DATA/$2 | cut -d\  -f2)		# Check if the second argument is txt file contating list of subjects
		if [[ $isFile = "ASCII" ]];then
			subjects=$(cat $DATA/$2)
		else
			subjects=$2
		fi

		for SUB in $subjects;do
			if [ ${SUB: -1} == "/" ];then
				SUB=${SUB%/}			# delete slash at the end of SUB varaible if is containted
			fi

			check_input d ${DATA}/${SUB}

			main $DATA $SUB

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
	show "\tData directory: $DATA" y
	show "\tSubject: $SUB" y
	show "\tMachine name: $HOSTNAME" y
	show "\tNum of CPU: $(nproc)" y
	show "\tSCT version: $SCT_VER" y
	show "$LINE" yt


	check_input b matlab fslmerge sct_dmri_moco

#-----------------------------------------------
# Diffusion preprocessing

	check_input d ${DATA}/${SUB}/Diffusion

	cd ${DATA}/${SUB}/Diffusion

	ls *nii.gz &>/dev/null		# check if Diff files are gunzipped
	if [[ $? = 2 ]];then
		show "Compressing diffusion data..."
		gzip *.nii
	fi

	# ZOOMit AP (either interp or notinterp)
	ZOOM_AP=`ls ${DATA}/${SUB}/Diffusion/*d_ep2d_diff_tra_3mm_ZOOMit_AP_63d_2sh_GN_DFC*.nii.gz | sed 's/.\{7\}$//'`
	for TEMP in $ZOOM_AP;do
		dimension=`${FSLDIR}/bin/fslval $TEMP pixdim1`
		if [ $dimension == 0.650000 ];then
			ZOOM_AP_INTERP=$TEMP
		elif [ $dimension == 1.300000 ];then
			ZOOM_AP_NOTINTERP=$TEMP
		fi
	done

	# ZOOMit PA (either interp or notinterp)
	ZOOM_PA=`ls ${DATA}/${SUB}/Diffusion/*d_ep2d_diff_tra_3mm_ZOOMit_PA_63d_2sh_GN_DFC*.nii.gz | sed 's/.\{7\}$//'`
	for TEMP in $ZOOM_PA;do
		dimension=`${FSLDIR}/bin/fslval $TEMP pixdim1`
		if [ $dimension == 0.650000 ];then
			ZOOM_PA_INTERP=$TEMP
		elif [ $dimension == 1.300000 ];then
			ZOOM_PA_NOTINTERP=$TEMP
		fi
	done


	#READOUT=`matlab -nosplash -nodisplay -nodesktop -nojvm -r "get_readout('_2007B_04d_ep2d_diff_tra_3mm_ZOOMit_AP_63d_2sh_GN_DFC_MIX_dicom_header.mat'),exit" | tail -c 20`	# call matlab function for calculation of Total Read Out time; I have computed times yet, So you can use constants belows


	for FOLDER in $ACQ;do
		if [ ! -d ${DATA}/${SUB}/Results/${FOLDER} ]; then
			#if [ `echo $FOLDER | tail -c 5` == moco ];then
			if [ $FOLDER == Diff_preproc_ZOOMit_interp ] || [ $FOLDER == Diff_preproc_ZOOMit_interp_moco ] ;then
				READOUT=0.0484
				diff_prep $ZOOM_AP_INTERP $ZOOM_PA_INTERP $READOUT $FOLDER	# Call diff_preop function
			elif [ $FOLDER == Diff_preproc_ZOOMit_notinterp ] || [ $FOLDER == Diff_preproc_ZOOMit_notinterp_moco ] ;then
				READOUT=0.0967
				diff_prep $ZOOM_AP_NOTINTERP $ZOOM_PA_NOTINTERP $READOUT $FOLDER	# Call diff_preop function
			fi
		else
			show "Preprocessing of diffusion data in ${DATA}/${SUB}/Results/${FOLDER} folder has been done before." g
		fi
	done


	if [[ $RESOLVE_RUN == 1 ]];then
		# RESOLVE (without interpolation) (seq 08 and 11)
		RESOLVE_AP=`ls ${DATA}/${SUB}/Diffusion/*d_resolve_AX_AP_p2.nii.gz | sed 's/.\{7\}$//'`
		RESOLVE_PA=`ls ${DATA}/${SUB}/Diffusion/*d_resolve_AX_PA_p2.nii.gz | sed 's/.\{7\}$//'`
		READOUT=0.0198
		FOLDER="Diff_preproc_RESOLVE"
		if [ ! -d ${DATA}/${SUB}/Results/${FOLDER} ]; then
			diff_prep_resolve $RESOLVE_AP $RESOLVE_PA $READOUT $FOLDER	# Call diff_preop_resolve function
		else
			show "Preprocessing of diffusion data in ${DATA}/${SUB}/Results/${FOLDER} folder has been done before." g
		fi
	fi

#-----------------------------------------------
# topup
	# ZOOMit
	NAME="zoomit"

	for FOLDER in $ACQ;do
		if [ ! -f ${DATA}/${SUB}/Results/${FOLDER}/b0_${NAME}_topup.nii.gz ]; then
			topup_function $DATA $FOLDER $NAME &		# Call function for topup
			pidlist="$pidlist $!"    # add process ID of most recent executed background pipeline
		else
			show "topup on data in ${DATA}/${SUB}/Results/${FOLDER} folder has been done before." g
		fi
	done

	if [[ $RESOLVE_RUN == 1 ]];then
		# RESOLVE (seq 08 and 11)
		FOLDER="Diff_preproc_RESOLVE"
		NAME="resolve"

		if [ ! -f ${DATA}/${SUB}/Results/${FOLDER}/b0_${NAME}_topup.nii.gz ]; then
			topup_function $DATA $FOLDER $NAME &		# Call function for topup
			pidlist="$pidlist $!"    # add process ID of most recent executed background pipeline
		else
			show "topup on data in ${DATA}/${SUB}/Results/${FOLDER} folder has been done before." g
		fi
	fi

	if [[ "$pidlist" != "" ]]; then
		wait $pidlist
		pidlist=""
	fi


#-----------------------------------------------
# eddy
	#ZOOMit
	NAME="zoomit"

	for FOLDER in $ACQ;do
		if [ ! -f ${DATA}/${SUB}/Results/${FOLDER}/eddy_${NAME}.nii.gz ]; then
			eddy_function $DATA $FOLDER $NAME &		# Call function for eddy
			pidlist="$pidlist $!"    # add process ID of most recent executed background pipeline
		else
			show "eddy in $FOLDER folder has been done before." g
		fi
	done

	if [[ $RESOLVE_RUN == 1 ]];then
		#RESOLVE
		NAME="resolve"
		FOLDER="Diff_preproc_RESOLVE"
		if [ ! -f ${DATA}/${SUB}/Results/${FOLDER}/eddy_${NAME}.nii.gz ]; then
			eddy_function $DATA $FOLDER $NAME &		# Call function for eddy
			pidlist="$pidlist $!"    # add process ID of most recent executed background pipeline
		else
			show "eddy on data in ${DATA}/${SUB}/Results/${FOLDER} folder has been done before." g
		fi
	fi

	if [[ "$pidlist" != "" ]]; then
		wait $pidlist
		pidlist=""
	fi


#----------------------------------------------
# ditfit
	#ZOOMit
	NAME="zoomit"

	for FOLDER in $ACQ;do
		if [ ! -f ${DATA}/${SUB}/Results/${FOLDER}/dMRI/dti_${NAME}_FA.nii.gz ]; then
			dtifit_function $DATA $FOLDER $NAME &		# Call function for dtifit
			pidlist="$pidlist $!"    # add process ID of most recent executed background pipeline
		else
			show "Estimation of DTI model using dtifit on data in $FOLDER folder is done." g
		fi
	done

	if [[ $RESOLVE_RUN == 1 ]];then
		#RESOLVE
		NAME="resolve"
		FOLDER="Diff_preproc_RESOLVE"
		if [ ! -f ${DATA}/${SUB}/Results/${FOLDER}/dMRI/dti_${NAME}_FA.nii.gz ]; then
			dtifit_function $DATA $FOLDER $NAME &		# Call function for dtifit
			pidlist="$pidlist $!"    # add process ID of most recent executed background pipeline
		else
			show "Estimation of DTI model using dtifit on data in $FOLDER folder is done." g
		fi
	fi

	if [[ "$pidlist" != "" ]]; then
		wait $pidlist
		pidlist=""
	fi
#----------------------------------------------
# registration between T2TRA and DIFF spaces (T2TRA to DIFF)
	#ZOOMit
	NAME="zoomit"

	for FOLDER in $ACQ;do
		if [ ! -f ${DATA}/${SUB}/Results/${FOLDER}/T2TRA_in_DIFF.nii.gz ]; then
			t2tra_to_diff_function $DATA $FOLDER $NAME &		# Call function for registration between T2TRA and DIFF spaces
			pidlist="$pidlist $!"    # add process ID of most recent executed background pipeline
		else
			show "Registration between T2TRA and DIFF spaces in $FOLDER folder has been done before." g
		fi
	done

	if [[ $RESOLVE_RUN == 1 ]];then
		#RESOLVE
		NAME="resolve"
		FOLDER="Diff_preproc_RESOLVE"
		if [ ! -f ${DATA}/${SUB}/Results/${FOLDER}/T2TRA_in_DIFF.nii.gz ]; then
			t2tra_to_diff_function $DATA $FOLDER $NAME &		# Call function for registration between T2TRA and DIFF spaces
			pidlist="$pidlist $!"    # add process ID of most recent executed background pipeline
		else
			show "Registration between T2TRA and DIFF spaces in $FOLDER folder has been done before." g
		fi
	fi

	if [[ "$pidlist" != "" ]]; then
		wait $pidlist
		pidlist=""
	fi
#----------------------------------------------
# diff labeling
	#ZOOMit
	NAME="zoomit"

	for FOLDER in $ACQ;do
		if [ ! -f ${DATA}/${SUB}/Results/${FOLDER}/dwi_mean_seg_labeled.nii.gz ]; then
			difflabel_function $DATA $FOLDER $NAME &		# Call function for labeling of diff data
			pidlist="$pidlist $!"    # add process ID of most recent executed background pipeline
	else
			show "Labeling of diff data in $FOLDER folder has been done before." g
		fi
	done

	if [[ $RESOLVE_RUN == 1 ]];then
		#RESOLVE
		NAME="resolve"
		FOLDER="Diff_preproc_RESOLVE"
		if [ ! -f ${DATA}/${SUB}/Results/${FOLDER}/dwi_mean_seg_labeled.nii.gz ]; then
			difflabel_function $DATA $FOLDER $NAME &		# Call function for labeling of diff data
			pidlist="$pidlist $!"    # add process ID of most recent executed background pipeline
		else
			show "Labeling of diff data in $FOLDER folder has been done before." g
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
		if [ ! -f ${DATA}/${SUB}/Results/${FOLDER}/dMRI/dti_${NAME}_FA_cord.nii* ]; then
			masking $DATA $FOLDER $NAME &		# Call function for masking
			pidlist="$pidlist $!"    # add process ID of most recent executed background pipeline
		else
			show "Masking off diffusion data and warping diffusion data into T2TRA in $FOLDER folder has been done before." g
		fi
	done

	if [[ $RESOLVE_RUN == 1 ]];then
		#RESOLVE
		NAME="resolve"
		FOLDER="Diff_preproc_RESOLVE"
		if [ ! -f ${DATA}/${SUB}/Results/${FOLDER}/dMRI/dti_${NAME}_FA_cord.nii* ]; then
			masking $DATA $FOLDER $NAME &		# Call function for masking
			pidlist="$pidlist $!"    # add process ID of most recent executed background pipeline
		else
			show "Masking off diffusion data and warping diffusion data into T2TRA in $FOLDER folder has been done before." g
		fi
	fi

	if [[ "$pidlist" != "" ]]; then
		wait $pidlist
		pidlist=""
	fi

#----------------------------------------------
# Permissions change

	permission=$(ls -ld ${DATA}/${SUB}/Results)
	if [ ${permission:1:3} != ${permission:4:3} ];then
		exe "chmod -R g=u ${DATA}/${SUB}/Results" v
	fi
	if [[ ${permission:7:3} =~ "rwx" ]];then
		exe "chmod -R o-rwx ${DATA}/${SUB}/Results" v
	fi

#----------------------------------------------
	show "$LINE" yt
	show "\t$SCRIPT_NAME for ${DATA}/${SUB} finished." y
	show "$LINE" yt

	# END of main function

}


#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Function which is called for diffusion preprocessing of ZOOMit data
diff_prep()
{
	ZOOM_AP=$1
	ZOOM_PA=$2
	READOUT=$3
	FOLDER=$4

	show "Starting preprocessing of diffusion data in ${DATA}/${SUB}/Results/${FOLDER} folder!" y

	check_input dc ${DATA}/${SUB}/Results			# Create Results folder
	check_input dc ${DATA}/${SUB}/Results/${FOLDER}		# Create Diff_preproc_ZOOMit_<> folder

	cd ${DATA}/${SUB}/Results/${FOLDER}

	# Create bval and bvec file for PA acquision (this acquision doesnt contain these files)
	bvalues=`fslval ${ZOOM_PA}.nii.gz dim4`				# Get number of b0 images in PA data
	indx="0"
	for ((i=1; i<$bvalues; ++i)); do indx="$indx 0"; done		# Create variable indx containing same number of zeros as number of b0 images in PA data
	echo $indx > temp.bval						# Create temp.bval file
	echo $indx > temp.bvec; echo $indx >> temp.bvec; echo $indx >> temp.bvec	# Create temp.bvec file (three lines file)

	echo `cat ${ZOOM_AP}.bval``cat temp.bval` > eddy_input_zoomit.bval	# Create eddy_input_zoomit.bval file (add temp.bval containing zeros to end of ${ZOOM_AP}.bval)
	paste -d "" ${ZOOM_AP}.bvec temp.bvec > eddy_input_zoomit.bvec		# Create eddy_input_zoomit.bvec file (add temp.bvec containing zeros to ends of ${ZOOM_AP}.bvec)
	rm temp.*

	exe "fslmerge -a eddy_input_zoomit.nii.gz ${ZOOM_AP}.nii.gz ${ZOOM_PA}.nii.gz" v	# Merge AP and PA ZOOMit data into one file
	if [ `echo $FOLDER | tail -c 5` == moco ];then
		show "Starting sct_dmri_moco in $FOLDER folder." y

		exe "sct_maths -i eddy_input_zoomit.nii.gz -mean t -o dmri_mean.nii.gz" v
		exe "sct_register_multimodal -i ${DATA}/${SUB}/Results/Anat_Preproc/T2TRA_thr_bias_corr_seg.nii.gz -d dmri_mean.nii.gz -identity 1 -x nn" v
		exe "sct_create_mask -i dmri_mean.nii.gz -p centerline,T2TRA_thr_bias_corr_seg_reg.nii.gz -size 110mm" v
		exe "sct_crop_image -i eddy_input_zoomit.nii.gz -m mask_dmri_mean.nii.gz -o dmri_crop.nii.gz" v
		exe "sct_dmri_moco -i dmri_crop.nii.gz -bvec eddy_input_zoomit.bvec" v

		exe "mv eddy_input_zoomit.nii.gz eddy_input_zoomit_without_moco.nii.gz" v
		exe "mv dmri_crop_moco.nii.gz ./eddy_input_zoomit.nii.gz" v
		#gzip eddy_input_zoomit.nii
		#rm -r tmp*/
		show "sct_dmri_moco in $FOLDER folder is done." y
	fi

	B0INDEX=""
	fslsplit eddy_input_zoomit.nii.gz temporary	# Split ZOOMit data into infividual images
	B0IND=0
	IMAGEID=0

	MERGECOMMAND="fslmerge -a b0_zoomit.nii.gz"
	for BVAL in `cat eddy_input_zoomit.bval`;do			# Create variable B0INDEX containg order of b0 images in merged file
		if [ $BVAL -eq 0 ] || [ $BVAL -eq 5 ];then
			B0IND=$(($B0IND+1))
			MERGECOMMAND="$MERGECOMMAND temporary`printf %04d $IMAGEID`.nii.gz"
		fi
		B0INDEX="$B0INDEX $B0IND"
		IMAGEID=$(($IMAGEID+1))
	done
	echo $B0INDEX > index.txt
	B0NUM=`cat index.txt | awk '{print $NF}'`			# Count number of b0 images in merged file
	B0AP=$(($B0NUM-$bvalues))					# Get number of b0 images in AP data

	echo "0 -1 0 $READOUT" > acq_file_zoomit.txt
	for MES in `seq 2 $B0NUM`;do					# Create acq_file which is necessary for topup
		if [ $MES -le $B0AP ];then
			echo "0 -1 0 $READOUT" >> acq_file_zoomit.txt
		else
			echo "0 1 0 $READOUT" >> acq_file_zoomit.txt
		fi
	done

	`echo $MERGECOMMAND`						# Create file containing only b0 images
	rm temporary*.nii.gz

	show "Preprocessing of diffusion data in ${DATA}/${SUB}/Results/${FOLDER} folder is done." y
}
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Function which is called for diffusion preprocessing of RESOLVE data
diff_prep_resolve()
{

	RESOLVE_AP=$1
	RESOLVE_PA=$2
	READOUT=$3
	FOLDER=$4

	show "Starting preprocessing of diffusion data in ${DATA}/${SUB}/Results/${FOLDER} folder!" y

	check_input dc ${DATA}/${SUB}/Results			# Create Results folder
	check_input dc ${DATA}/${SUB}/Results/${FOLDER}		# Create Diff_preproc_RESOLVE folder

	cd ${DATA}/${SUB}/Results/${FOLDER}

	exe "fslmerge -a eddy_input_resolve.nii.gz ${RESOLVE_AP}.nii.gz ${RESOLVE_PA}.nii.gz" v	# Merge AP and PA ZOOMit data into one file

	echo `cat ${RESOLVE_AP}.bval``cat ${RESOLVE_PA}.bval` > eddy_input_resolve.bval	# Create eddy_input_resolve.bval file by merge bval files of RESOLVE_AP and RESOLVE_PA
	paste -d "" ${RESOLVE_AP}.bvec ${RESOLVE_PA}.bvec > eddy_input_resolve.bvec	# Create eddy_input_resolve.bvec file by merge bvec files of RESOLVE_AP and RESOLVE_PA


	B0INDEX=""
	fslsplit eddy_input_resolve.nii.gz temporary	# Split RESOLVE data into infividual images
	B0IND=0
	IMAGEID=0

	MERGECOMMAND="fslmerge -a b0_resolve.nii.gz"
	for BVAL in `cat eddy_input_resolve.bval`;do			# Create variable B0INDEX containg order of b0 images in merged file
		if [ $BVAL -eq 0 ];then
			B0IND=$(($B0IND+1))
			MERGECOMMAND="$MERGECOMMAND temporary`printf %04d $IMAGEID`.nii.gz"
		fi
		B0INDEX="$B0INDEX $B0IND"
		IMAGEID=$(($IMAGEID+1))
	done
	echo $B0INDEX > index.txt
	B0NUM=`cat index.txt | awk '{print $NF}'`			# Count number of b0 images in merged file
	B0AP=$(($B0NUM/2))

	echo "0 -1 0 $READOUT" > acq_file_resolve.txt
	for MES in `seq 2 $B0NUM`;do					# Create acq_file which is necessary for topup
		if [ $MES -le $B0AP ];then
			echo "0 -1 0 $READOUT" >> acq_file_resolve.txt
		else
			echo "0 1 0 $READOUT" >> acq_file_resolve.txt
		fi
	done

	`echo $MERGECOMMAND`						# Create file containing only b0 images
	rm temporary*.nii.gz

	show "Preprocessing of diffusion data in ${DATA}/${SUB}/Results/${FOLDER} folder is done." y

}
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Function which is called for topup and mean of data after topup
topup_function()
{
	DATA=$1
	FOLDER=$2
	NAME=$3

	show "Starting topup on data in $FOLDER folder!" y

	cd ${DATA}/${SUB}/Results/${FOLDER}

	exe "topup --imain=b0_${NAME}.nii.gz --datain=acq_file_${NAME}.txt --out=topup_${NAME} --subsamp=1,1,1,1,1,1,1,1,1 --config=b02b0.cnf --iout=b0_${NAME}_topup --fout=field_${NAME}_topup" v

	exe "fslmaths b0_${NAME}_topup.nii.gz -Tmean b0_${NAME}_topup_mean.nii.gz" v	# mean b0

	show "topup on data in ${DATA}/${SUB}/Results/${FOLDER} folder is done." y

}
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Function which is called for eddy (eddy uses either manually segmented mask of SC or whole binarized FOV)
eddy_function()
{
	DATA=$1
	FOLDER=$2
	NAME=$3

	show "Starting eddy on data in $FOLDER folder!" y

	cd ${DATA}/${SUB}/Results/${FOLDER}

	# Note: originally the binary was called eddy
	if [[ "$(which eddy)" != "" ]]; then
		EDDYCMD="eddy"
	elif [[ "$(which eddy_openmp)" != "" ]]; then
		EDDYCMD="eddy_openmp"
	else
		show "eddy binary is not found or path is wrong. Fix it!" e
	fi


	if [ -f ${DATA}/${SUB}/Results/${FOLDER}/b0_${NAME}_topup_mean_seg.nii.gz ];then
		show "Starting eddy on data in ${DATA}/${SUB}/Results/${FOLDER} folder with manually created mask!" y
		exe "${EDDYCMD} --imain=eddy_input_${NAME}.nii.gz --mask=b0_${NAME}_topup_mean_seg.nii.gz --index=index.txt --acqp=acq_file_${NAME}.txt --bvecs=eddy_input_${NAME}.bvec --bvals=eddy_input_${NAME}.bval --topup=topup_${NAME} --out=eddy_${NAME} -v" v
		show "eddy on data in ${DATA}/${SUB}/Results/${FOLDER} folder with manually created mask is done!" y
	elif [ ! -f ${DATA}/${SUB}/Results/${FOLDER}/b0_${NAME}_topup_mean_seg.nii.gz ] && [ ! -f ${DATA}/${SUB}/Results/${FOLDER}/b0_${NAME}_topup_mean_bin.nii.gz ];then
		show "Manually segmented mask of the spinal cord in ${DATA}/${SUB}/Results/${FOLDER} with name b0_${NAME}_topup_mean_seg.nii.gz is not found!" y
		show "Starting eddy on data in ${DATA}/${SUB}/Results/${FOLDER} folder with binarized FOV as mask!" y
		exe "fslmaths b0_${NAME}_topup_mean.nii.gz -bin b0_${NAME}_topup_mean_bin.nii.gz" v
		exe "${EDDYCMD} --imain=eddy_input_${NAME}.nii.gz --mask=b0_${NAME}_topup_mean_bin.nii.gz --index=index.txt --acqp=acq_file_${NAME}.txt --bvecs=eddy_input_${NAME}.bvec --bvals=eddy_input_${NAME}.bval --topup=topup_${NAME} --out=eddy_${NAME} -v" v
		show "eddy on data in ${DATA}/${SUB}/Results/${FOLDER} folder with binarized FOV as mask is done!" y
	fi

}

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Function which is called for estimation of DTI model using dtifit
dtifit_function()
{
	DATA=$1
	FOLDER=$2
	NAME=$3

	show "Starting estimation of DTI model using dtifit on data in $FOLDER folder." y

	check_input dc ${DATA}/${SUB}/Results/${FOLDER}/dMRI
	cd ${DATA}/${SUB}/Results/${FOLDER}/

	exe "cp eddy_${NAME}.nii.gz dMRI/data.nii.gz" v
	exe "cp eddy_input_${NAME}.bval dMRI/bvals" v
	exe "cp eddy_${NAME}.eddy_rotated_bvecs dMRI/bvecs" v

	if [ -f ${DATA}/${SUB}/Results/${FOLDER}/b0_${NAME}_topup_mean_seg.nii.gz ];then
		exe "cp b0_${NAME}_topup_mean_seg.nii.gz dMRI/nodif_brain_mask.nii.gz" v
	elif [ ! -f ${DATA}/${SUB}/Results/${FOLDER}/b0_${NAME}_topup_mean_seg.nii.gz ] && [ -f ${DATA}/${SUB}/Results/${FOLDER}/b0_${NAME}_topup_mean_bin.nii.gz ];then
		exe "cp b0_${NAME}_topup_mean_bin.nii.gz dMRI/nodif_brain_mask.nii.gz" v
	else
		error "Problem during copy data into dMRI folder happened (b0_${NAME}_topup_mean_seg.nii.gz nor b0_${NAME}_topup_mean_bin.nii.gz dMRI/nodif_brain_mask.nii.gz do not exist)!"
	fi

	cd ${DATA}/${SUB}/Results/${FOLDER}/dMRI
	exe "dtifit -k data.nii.gz -o dti_${NAME} -m nodif_brain_mask.nii.gz -r bvecs -b bvals -w" v

	show "Estimation of DTI model using dtifit on data in $FOLDER folder is done." y

}

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Function which is called for T2TRA coregistration to DIFF space using sct_register_multimodal function
# Registration is done in two steps:
#		1) T2-w axial SC segmentation (T2TRA_thr_bias_corr_seg.nii.gz) is registred to DWI space
# 	2) Mean b0 image (b0_mean) is registred to T2-w axial (T2TRA) image. Segmentations of SC from both images are used for improving of final registration. Warping field from previous step is used.
t2tra_to_diff_function()
{
	DATA=$1
	FOLDER=$2
	NAME=$3

	show "Starting T2TRA and diffusion data coregistration in ${DATA}/${SUB}/Results/${FOLDER} folder." y

	cd ${DATA}/${SUB}/Results/${FOLDER}

	exe "cp eddy_${NAME}.eddy_rotated_bvecs bvecs" v			# necessary, because sct_dmri_separate_b0_and_dwi does not recognize file type .eddy_rotated_bvecs
	exe "sct_dmri_separate_b0_and_dwi -i eddy_${NAME}.nii.gz -bvec bvecs -a 1 -bval eddy_input_${NAME}.bval -r 1" v
	exe "rm bvecs" v
	exe "sct_register_multimodal -i ${DATA}/${SUB}/Results/Anat_Preproc/T2TRA_thr_bias_corr_seg.nii.gz -d dwi_mean.nii.gz -identity 1 -x nn -o T2TRA_seg_in_DIFF.nii.gz" v
	exe "mv warp_dwi_mean2T2TRA_thr_bias_corr_seg.nii.gz warp_DIFFtoT2TRA.nii.gz" v
	exe "mv warp_T2TRA_thr_bias_corr_seg2dwi_mean.nii.gz warp_T2TRAtoDIFF.nii.gz" v
	exe "mv T2TRA_seg_in_DIFF_inv.nii.gz dwi_mean_in_T2TRA.nii.gz" v
	if [ $SUB == 2180B ];then							# SUB 2180B has so large scale of intensities in dwi_mean image -> necessary to threshold it, JV, 7.5.18
		exe "fslmaths dwi_mean.nii.gz -uthr 500 dwi_mean_thr.nii.gz" v
		exe "sct_propseg -i dwi_mean_thr.nii.gz -c dwi -init-centerline T2TRA_seg_in_DIFF.nii.gz -radius 6 -max-deformation 5" v
		exe "mv dwi_mean_thr_seg.nii.gz dwi_mean_seg.nii.gz" v
		exe "mv dwi_mean_thr_centerline.nii.gz dwi_mean_centerline.nii.gz" v
	else

		if [ $DEEPSEG == 0 ];then

			show "Starting SC segmentation in ${DATA}/${SUB}/Results/${FOLDER} folder using sct_propseg." y
			exe "sct_propseg -i dwi_mean.nii.gz -c dwi -init-centerline T2TRA_seg_in_DIFF.nii.gz -radius 6 -max-deformation 5" v
			show "SC segmentation in ${DATA}/${SUB}/Results/${FOLDER} folder using sct_propseg is done." y

		elif [ $DEEPSEG == 1 ];then

			show "Starting SC segmentation in ${DATA}/${SUB}/Results/${FOLDER} folder using sct_deepseg." y
			exe "sct_deepseg_sc -i dwi_mean.nii.gz -c dwi" v
			show "SC segmentation in ${DATA}/${SUB}/Results/${FOLDER} folder using sct_deepseg is done." y

		fi
	fi
	# Register DIFF to T2TRA
	exe "sct_register_multimodal -i b0_mean.nii.gz -d ${DATA}/${SUB}/Results/Anat_Preproc/T2TRA_thr_bias_corr.nii.gz -iseg dwi_mean_seg.nii.gz -dseg ${DATA}/${SUB}/Results/Anat_Preproc/T2TRA_thr_bias_corr_seg.nii.gz -param step=1,type=seg,algo=slicereg,metric=CC,smooth=5:step=2,type=seg,algo=bsplinesyn,metric=CC,smooth=1,iter=3 -initwarp warp_DIFFtoT2TRA.nii.gz -initwarpinv warp_T2TRAtoDIFF.nii.gz" v
	exe "mv T2TRA_thr_bias_corr_reg.nii.gz T2TRA_in_DIFF.nii.gz" v
	exe "mv b0_mean_reg.nii.gz b0_mean_in_T2TRA.nii.gz" v
	exe "mv warp_T2TRA_thr_bias_corr2b0_mean.nii.gz warp_T2TRAtoDIFF.nii.gz" v
	exe "mv warp_b0_mean2T2TRA_thr_bias_corr.nii.gz warp_DIFFtoT2TRA.nii.gz" v
	#exe "rm -r qc" v

	show "T2TRA and diffusion data coregistration in $FOLDER folder is done!" y
}
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Masking of results from dtifit and warping these results into T2TRA space
masking()
{
	DATA=$1
	FOLDER=$2
	NAME=$3

	show "Starting masking of diffusion data in ${DATA}/${SUB}/Results/${FOLDER} folder." y

	cd ${DATA}/${SUB}/Results/${FOLDER}

	exe "fslmaths dMRI/dti_${NAME}_FA.nii.gz -mas dwi_mean_seg.nii.gz dMRI/dti_${NAME}_FA_cord.nii.gz" v
	#gunzip dMRI/dti_${NAME}_FA_cord.nii.gz
	exe "fslmaths T2TRA_in_DIFF.nii.gz -mas dwi_mean_seg.nii.gz T2TRA_in_DIFF_cord.nii.gz" v
	#gunzip T2TRA_in_DIFF_cord.nii.gz
	exe "fslmaths dMRI/dti_${NAME}_MD.nii.gz -mas dwi_mean_seg.nii.gz dMRI/dti_${NAME}_MD_cord.nii.gz" v
	#gunzip dMRI/dti_${NAME}_MD_cord.nii.gz

	show "masking of diffusion data in $FOLDER folder is done!" y

# Warping DIFF results into T2TRA
	show "Starting warping diffusion data into T2TRA in ${DATA}/${SUB}/Results/${FOLDER} folder." y

	check_input dc ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra

	for j in dMRI/dti_${NAME}_FA.nii.gz dMRI/dti_${NAME}_FA_cord.nii.gz dMRI/dti_${NAME}_MD.nii.gz dMRI/dti_${NAME}_MD_cord.nii.gz dwi_mean_seg.nii.gz b0_mean.nii.gz dwi_mean.nii.gz dMRI/dti_${NAME}_V1.nii.gz dMRI/dti_${NAME}_L1.nii.gz dMRI/dti_${NAME}_L2.nii.gz dMRI/dti_${NAME}_L3.nii.gz field_${NAME}_topup.nii.gz;do
		OUTPUT=`echo $j | sed 's:.*/::g'`	#delete string before slash if is contained
		exe "sct_apply_transfo -i ${DATA}/${SUB}/Results/${FOLDER}/$j -d ${DATA}/${SUB}/Results/Anat_Preproc/T2TRA_thr_bias_corr.nii.gz -w ${DATA}/${SUB}/Results/${FOLDER}/warp_DIFFtoT2TRA.nii.gz -o ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra/$OUTPUT" v
	done

	exe "fslmaths ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra/dwi_mean_seg.nii.gz -thr 0.75 -bin ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra/dwi_mean_seg.nii.gz" v
	exe "fslmaths ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra/dti_${NAME}_FA_cord.nii.gz -thr 0.1 ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra/dti_${NAME}_FA_cord.nii.gz" v
	exe "fslmaths ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra/dti_${NAME}_MD_cord.nii.gz -thr 0.00001 -mas ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra/dwi_mean_seg.nii.gz ${DATA}/${SUB}/Results/${FOLDER}/diff_in_t2tra/dti_${NAME}_MD_cord.nii.gz" v
}
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Function which is called for vertebrae levels labeling in diff space
difflabel_function()
{
	DATA=$1
	FOLDER=$2
	NAME=$3

	show "Starting labeling of diffusion data in ${DATA}/${SUB}/Results/${FOLDER} folder!" y

	cd ${DATA}/${SUB}/Results/${FOLDER}

	exe "sct_apply_transfo -i ${DATA}/${SUB}/Results/Anat_Preproc/T2TRA_thr_bias_corr_seg_labeled.nii.gz -d b0_mean.nii.gz -w warp_T2TRAtoDIFF.nii.gz -o b0_mean_warped_labels.nii.gz -x nn" v
	exe "fslmaths b0_mean_warped_labels.nii.gz -mas dwi_mean_centerline.nii.gz -Xmax -Ymax -uthr 3 -thr 3 output" v
	C3C4SLICE=`fslstats output.nii.gz -w | cut -d' ' -f5`
	exe "rm output.nii.gz" v
	exe "sct_label_vertebrae -i b0_mean.nii.gz -s dwi_mean_seg.ni* -c t2 -initz $C3C4SLICE,3 -r 1" v
	#sct_label_vertebrae -i T2TRA_in_DIFF.nii.gz -s dwi_mean_seg.ni* -c t2 -initz $C3C4SLICE,3 -r 1

	show "Labeling of diffusion data in $FOLDER folder is done!" y
}

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

init $@
