#!/bin/bash

#-----------------------------------------------------------------------------------------------------------
# Script for processing of the spinal cord data
# USAGE:
# 		run_analysis_final.sh <DATA_FOLDER> <SUB_ID> <SEQ_order>
#						- DATA_FOLDER - folder containing individual subjects data
#						- SUB_ID - name of the subject folder with raw (nifti) anatomical and diffusion data
#						- SEQORDER - order of diffusion protocols (ZOOMit_interp, ZOOMit_interp_moco, ZOOMit_nointerp, ZOOMit_nointerp_moco, RESOLVE)
# EXAMPLE:
#			run_analysis_final.sh /home/user/data 1234B 11001
#
# Analysis contains:
#		1) Bias field correction of orignal T2TRA and T2SAG images using N4BiasFieldCorrection tool (part of ANTS) and thresholding of low values
# 	2) Slice-by-slice correction of zig-zag artifact of T2TRA image (could be switched off inside script)
# 	3) Resampling and cropping of T2 sagittal image (T2SAG) (for fitting with resolution of T2TRA image)
# 	4) Spinal cord segmentation of T2 sagittal image (T2SAG) either by sct_deepseg_sc or sct_propseg function (could be set inside script)
# 	5) Vertebrae labeling of T2SAG
# 	6) Register T2SAG to template (T2SAG SC segmentation is necessary)
# 	7) Spinal cord segmentation of T2TRA etiher by sct_propseg or using T2SAGseg as init-centerline or manually created init-centerline) or by sct_deepseg_sc (could be set inside script)
# 	8) Segmentation of T2TRA GM either by sct_segment_graymatter or by sct_deepseg_gm
# 	9) Labeling of T2TRA SC and its correction at segment edges
# 	10) Register template and atlas (PAM50) to T2 transversal image (T2TRA space) using information from previous template to T2SAG registration
# 	11) Improve template registration to T2TRA by acounting WM/GM shape (requires correct WM and GM segmentation of T2TRA image)
#
# Jan Valosek, Rene Labounek. fMRI laboratory Olomouc, Czech Republic, 2017-2020
#
#-----------------------------------------------------------------------------------------------------------

# Some global variable
VERSION=25-02-2020

SCT_VER_GIT="7af1436782a5ea87758ac893a5f658d58eb3c60a"		# SCT v3.2.3 - https://github.com/neuropoly/spinalcordtoolbox/commit/7af1436782a5ea87758ac893a5f658d58eb3c60a

SCRIPT_NAME=${0##*/}													# get script name
SCRIPT_DIR="$(cd "$(dirname $0)";pwd -P)"			# get directory where is this script located
LINE="========================================================================================================================================="

DEEPSEG=1					# Type of SC segmentation: 0 - sct_propseg, 1 - sct_deepseg_sc
SLICE_THIC=2.5		# Slice thickness of T2TRA (2.5mm or 3mm)
RUNSLICE=1				# Slice-by-slice correction
BAD_T2SAG_SEG=0		# Result of T2SAG segmentation (if bad, the variable is set to 1 and T2SAG labeling, registration to template and T2 fusing is skipped)

NIFTI=${SCRIPT_DIR}/NIfTI_tools					# nifti tools for matlab
BINDIR=${SCRIPT_DIR}										# our bash scripts
SCTDIR=/usr/local/lib/sct								# path to SCT
ANTSDIR=/usr/local/lib/ants_v2.1.0			# BiasFiledCorrection path (part of ANTS)
PATH=$ANTSDIR:$BINDIR:$PATH

# Trap kill from check_input function if SIGUSR1 is set
trap "echo Exiting>&2;exit" SIGUSR1

print_help()
{

	echo -e "$LINE\nHelp for script performing anatomical analysis of the spinal cord data, version $VERSION."
	echo -e "USAGE:\n\t$SCRIPT_NAME <DATA_FOLDER> <SUB_ID> <SEQ_order>\nexample:\n\t$SCRIPT_NAME /md2/CSD 2256B 11001"
	echo -e "\nAnalysis contains folowing steps:"
	echo -e "\t1) Bias field correction of orignal T2TRA and T2SAG images using N4BiasFieldCorrection tool (part of ANTS) and thresholding of low values"
	echo -e "\t2) Slice-by-slice correction of zig-zag artifact of T2TRA image (could be switched off inside script)"
	echo -e "\t3) Resampling and cropping of T2 sagittal image (T2SAG) (for fitting with resolution of T2TRA image)"
	echo -e "\t4) Spinal cord segmentation of T2 sagittal image (T2SAG) either by sct_deepseg_sc or sct_propseg function (could be set inside script)"
	echo -e "\t5) Vertebrae labeling of T2SAG"
	echo -e "\t6) Register T2SAG to template (T2SAG SC segmentation is necessary)"
	echo -e "\t7) Spinal cord segmentation of T2TRA etiher by sct_propseg or using T2SAGseg as init-centerline or manually created init-centerline) or by sct_deepseg_sc (could be set inside script)"
	echo -e "\t8) Segmentation of T2TRA GM either by sct_segment_graymatter or by sct_deepseg_gm"
	echo -e "\t9) Labeling of T2TRA SC and its correction at segment edges"
	echo -e "\t10) Register template and atlas (PAM50) to T2 transversal image (T2TRA space) using information from previous template to T2SAG registration"
	echo -e "\t11) Improve template registration to T2TRA by acounting WM/GM shape (requires correct WM and GM segmentation of T2TRA image)"
	echo -e "\nValosek, Labounek 2017-2020\tfMRI laboratory, Olomouc, CZ\n$LINE"
	exit


}

# function for labeling of T2 transversal/axial (T2TRA) image
# labeling algorithm uses either information from T2SAG labeling or manually entered C2/3 slice
# if C2/3 slice is entered manually, it is saved into vert_index_tra.txt file
t2tra_labeling()
{
	# Perform labeling of T2TRA for the first time (as the initial information for C3/4 slice is used information from labeling of T2SAG image)
	if [[ `cat ${BINDIR}/vert_index_tra.txt | grep $SUB` = "" ]] && [[ $BAD_T2SAG_SEG = 0 ]];then

		exe "sct_apply_transfo -i T2SAG_resample_seg_labeled.nii.gz -d T2TRA_thr_bias_corr.nii.gz -w warp_T2SAG_resample_seg2T2TRA_thr_bias_corr.nii.gz -o T2TRA_warped_labels.nii.gz -x nn" v
		exe "fslmaths T2TRA_warped_labels.nii.gz -Xmax -Ymax -uthr 3 -thr 3 output" v
		C3C4SLICE=`fslstats output.nii.gz -w | cut -d' ' -f5`
		exe "rm output.nii.gz" v
		exe "sct_label_vertebrae -i T2TRA_thr_bias_corr.nii.gz -s T2TRA_thr_bias_corr_seg.nii.gz -c t2 -initz $C3C4SLICE,3 -r 1" v

	# Perform labeling of T2TRA for the first time with manually entered C2/3 slice (Segmentation of SC from T2SAG was not succesfull)
	elif [[ `cat ${BINDIR}/vert_index_tra.txt | grep $SUB` = "" ]] && [[ $BAD_T2SAG_SEG = 1 ]];then

		show "Labeling of T2TRA does not exist, open FSLeyes by following command:" w
		show "fsleyes ${DATA}/${SUB}/Results/Anat_Preproc/T2TRA_thr_bias_corr.nii.gz -dr 0 0.6 ${DATA}/${SUB}/Results/Anat_Preproc/T2TRA_thr_bias_corr_seg.nii.gz -a 50 -cm red &"
		show "Enter subjectID, number of C2/C3 slice and size of vertebrae in SI direction (eg. ${SUB} 45 10)): "
		return 1

	# Perform labeling based on information about C2/C3 level and vertebra size (in SI direction) from vert_index_tra.txt file
	else

		while read line; do
			if [[ $line = ${SUB}* ]];then
				LEVEL=$(cut -c 6-8 <<< $line);
				SIZE=$(cut -c 10- <<< $line);
				#if [[ $SIZE == "" ]];then SIZE=19;fi
			fi;
		done < ${BINDIR}/vert_index_tra.txt

		show "Computing labeling with manually entered parameters - C2/C3 level: $LEVEL and IS size: $SIZE"
		exe "sct_label_vertebrae -i T2TRA_thr_bias_corr.nii.gz -s T2TRA_thr_bias_corr_seg.nii.gz -c t2 -initz ${LEVEL},2 -param size_IS=${SIZE}" v

	fi

	show "Labeling and GM segmentation could be wrong - check it by FSLeyes by following command:" w
	show "fsleyes ${DATA}/${SUB}/Results/Anat_Preproc/T2TRA_thr_bias_corr.nii.gz -dr 0 0.6 ${DATA}/${SUB}/Results/Anat_Preproc/T2TRA_thr_bias_corr_gmseg.nii.gz -a 50 -cm red ${DATA}/${SUB}/Results/Anat_Preproc/T2TRA_thr_bias_corr_seg_labeled.nii.gz -a 50 -cm subcortical &" t
	show "Enter your answer and press [ENTER] (if labeling is ok: y/Y/a/A or if labeling is bad enter subjectID, number of C2/C3 slice and size of vertebrae in SI direction (eg. ${SUB} 45 10)): " t

}


init_task()
{

	# Print help
	if [ $# -eq 0 ] || [ $1 == "-h" ] || [ $# -ne 3 ];then
		print_help
	fi

# INICIALIZATION --------------------------------------------------------------------------
# Folders, Path etc.

	# Fetch some basic functions
	source ${BINDIR}/bash_basic_functions.sh

	DATA=$1
	if [ ${DATA: -1} == "/" ];then
		DATA=${DATA%/}			# delete slash at the end if is containted
	fi

	check_input d $DATA

	SUB=$2
	if [ ${SUB: -1} == "/" ];then
		SUB=${SUB%/}			# delete slash at the end of SUB varaible if is containted
	fi

	check_input d $DATA/$SUB

	SEQORDER=$3

	export LOGPATH=${DATA}/${SUB}/$(echo $SCRIPT_NAME | sed 's/.\{3\}$//').log
	exec > >(tee -a $LOGPATH) 2>&1

	# call function for checking requirments (MATLAB, FSL, N4BiasCorrection tool)
	check_input b matlab N4BiasFieldCorrection fslmerge

	if [[ "$(which sct_register_graymatter)" == "" ]];then
		if [[ "$(which sct_register_multimodal)" == "" ]];then
			show "SCT is not installed or path is wrong. Fix it!" e
		else
			show "SCT is installed but in version which does not contain sct_register_graymatter function (this function has been replaced in higger version of SCT by sct_register_multimodal)" e
		fi
	fi

	# Get SCT version
	SCT_VER=$(sct_check_dependencies | (grep "Spinal Cord Toolbox") | awk -F '/' '{print $NF}' | sed 's/)//')		# run sct_check_dependencies script, get line containing 'Spinal Cord Toolbox', delete everything before last slash and delete last bracket

	if [[ $SCT_VER != $SCT_VER_GIT ]]; then
		show "SCT v3.2.3 was not found. Current version of SCT is different - $SCT_VER\nDownload SCT v3.2.3 by:\ngit clone --branch v3.2.3 https://github.com/neuropoly/spinalcordtoolbox.git" e
	fi

	show "$LINE" yt
	show "\tStarting: $SCRIPT_NAME script, version: $VERSION" y
	show "\tData directory: $DATA" y
	show "\tSubject: $SUB" y
	show "\tDiff seq orded: $SEQORDER" y
	show "\tMachine name: $HOSTNAME" y
	show "\tNum of CPU: $(nproc)" y
	show "\tSCT version: $SCT_VER" y
	show "$LINE" yt

	if [ ! -d ${DATA}/${SUB}/Anat ];then
		show "Folder with anatomical data does not exist." e
	fi

	# Change permissions
	permission=$(ls -ld $DATA/$SUB/)
	if [ ${permission:1:3} != ${permission:4:3} ];then
		exe "chmod -R g=u $DATA/$SUB" v
	fi
	if [[ ${permission:7:3} =~ [rwx] ]];then
		exe "chmod -R o-rwx $DATA/$SUB" v
	fi

	cd ${DATA}/${SUB}/Anat

	# Check if input anatomical files (.nii) are compressed
	ls *nii.gz &>/dev/null
	if [[ $? = 2 ]];then
		show "Compressing anatomical data..."
		gzip *.nii
	fi

	# Fetch T2 sagital image
	T2SAG=`ls $DATA/$SUB/Anat/*s_t2_tse_sag_mm66_mm56_1mm3.nii.gz | sed 's/.\{7\}$//'` # delete last 7 characters (.nii.gz) from original file name

	# Fetch T2 axial/transversal image (could be 2.5 or 3 mm slice thickness)
	TEMP=`ls $DATA/$SUB/Anat/*s_t2_me2d_tra_p2.nii.gz | sed 's/.\{7\}$//'`

	for f in $TEMP;do
		dimension=$(fslval $f pixdim3)
		if [[ ${dimension:0:3} == $SLICE_THIC ]];then				# Select T2TRA with given slice-thin
			T2TRA=`ls $f.ni*`
		fi
	done

	if [ ! -d ${DATA}/${SUB}/Results ];then
		mkdir ${DATA}/${SUB}/Results
	fi

	if [ ! -d ${DATA}/${SUB}/Results/Anat_Preproc ]; then
		mkdir ${DATA}/${SUB}/Results/Anat_Preproc
		chmod 770 ${DATA}/${SUB}/Results/Anat_Preproc
	fi

	cd ${DATA}/${SUB}/Results/Anat_Preproc

	#-------------------------

	# call main function
	anat_analysis

	# function for analysis of diffusion data
 	diff_analysis


	#-------------------------

	# Permissions change
	permission=$(ls -ld $DATA/$SUB/Results)
	if [ ${permission:1:3} != ${permission:4:3} ];then
		exe "chmod -R g=u $DATA/$SUB/Results" v
	fi
	if [[ ${permission:7:3} =~ [rwx] ]];then
		exe "chmod -R o-rwx $DATA/$SUB/Results" v
	fi

	show "$LINE" yt
	show "\t$SCRIPT_NAME for $DATA/$SUB finished. Log is stored in $LOGPATH" y
	show "$LINE" yt

}


anat_analysis()
{
#----------------------------------------------------------------------------------------------------
# N4BiasFieldCorrection and thresholding of T2TRA and T2SAG (added by Jan Valosek 8.8.17)

	if [ ! -f T2SAG_thr_bias_corr.nii.gz ] || [ ! -f T2TRA_thr_bias_corr.nii.gz ];then

		show "Starting N4BiasFieldCorrection..." y

		# Tresholding
		if [ $SUB == 1860B ] || [ $SUB == 2684B ] || [ $SUB == 2825B ];then
			exe "fslmaths ${T2TRA}.nii.gz -thr 50 T2TRA_thr.nii.gz" v		#added for 1860B pacient; JV, 30.11.17
		elif [ $SUB == 2481B ];then
			exe "fslmaths ${T2TRA}.nii.gz -thr 30 T2TRA_thr.nii.gz" v		#added for 2481B ; JV, 01.10.18
		else
			exe "fslmaths ${T2TRA}.nii.gz -thr 70 T2TRA_thr.nii.gz" v
		fi
		exe "fslmaths ${T2SAG}.nii.gz -thr 20 T2SAG_thr.nii.gz" v
		# Binarize
		exe "fslmaths T2TRA_thr.nii.gz -bin T2TRA_thr_bin.nii.gz" v
		exe "fslmaths T2SAG_thr.nii.gz -bin T2SAG_thr_bin.nii.gz" v

		exe "gunzip -f T2TRA_thr_bin.nii.gz" v
		exe "gunzip -f T2SAG_thr_bin.nii.gz" v

		# Call m-file fill_holes.m for filling holes in thresholded binarized image and remove small noise areas
		exe "run_matlab addpath('$BINDIR','$NIFTI'),fill_holes('T2TRA_thr_bin'),exit" v
		exe "run_matlab addpath('$BINDIR','$NIFTI'),fill_holes('T2SAG_thr_bin'),exit" v

		exe "gzip -f T2TRA_thr_bin.nii T2TRA_thr_bin_fill.nii T2SAG_thr_bin.nii T2SAG_thr_bin_fill.nii" v

		# Multiplying thresholded binarized filled image with original image -> thresholded original without holes
		exe "fslmaths T2TRA_thr_bin_fill.nii.gz -mul $T2TRA.nii.gz T2TRA_thr_mul.nii.gz" v
		exe "fslmaths T2SAG_thr_bin_fill.nii.gz -mul $T2SAG.nii.gz T2SAG_thr_mul.nii.gz" v

		exe "gunzip -f T2TRA_thr_mul.nii.gz T2SAG_thr_mul.nii.gz" v

		# Bias correction
		exe "N4BiasFieldCorrection -d 3 -i T2TRA_thr_mul.nii -o T2TRA_thr_bias_corr.nii" v
		exe "N4BiasFieldCorrection -d 3 -i T2SAG_thr_mul.nii -o T2SAG_thr_bias_corr.nii" v

		exe "gzip -f T2TRA_thr_mul.nii T2SAG_thr_mul.nii T2TRA_thr_bias_corr.nii T2SAG_thr_bias_corr.nii" v

		# Normalization
		exe "fslmaths T2TRA_thr_bias_corr.nii.gz -div `fslstats T2TRA_thr_bias_corr.nii.gz -R | cut -f 2 -d " "` T2TRA_thr_bias_corr.nii.gz" v
		exe "fslmaths T2SAG_thr_bias_corr.nii.gz -div `fslstats T2SAG_thr_bias_corr.nii.gz -R | cut -f 2 -d " "` T2SAG_thr_bias_corr.nii.gz" v

		# Following images can be removed for save disk space
		exe "rm T2TRA_thr.nii.gz T2SAG_thr.nii.gz T2TRA_thr_bin.nii.gz T2SAG_thr_bin.nii.gz T2TRA_thr_bin_fill.nii.gz T2SAG_thr_bin_fill.nii.gz T2TRA_thr_mul.nii.gz T2SAG_thr_mul.nii.gz" v

		show "N4BiasFieldCorrection is done." y

		# Outputs for next use - T2SAG_thr_bias_corr.nii.gz and T2TRA_thr_bias_corr.nii.gz
	else
		show "N4BiasFieldCorrection has been done before." g
	fi

#----------------------------------------------------------------------------------------------------
# Slice-by-slice correction of zig-zag artifact of T2TRA image

	if [ $RUNSLICE == 1 ];then
		if [ ! -f isct_antsSliceRegularizedRegistration.log ];then

			show "Starting slice-by-slice correction of T2TRA image..." y

			cd ${DATA}/${SUB}/Results/Anat_Preproc

			exe "run_matlab addpath('$BINDIR','$NIFTI'),T2TRA_shift_correction_interp('T2TRA_thr_bias_corr','$BINDIR'),exit" v	# Call matlab script for slice-by-slice correction (based on split T2TRA into odd and even slices)

			exe "fslmaths T2TRA_thr_bias_corr_odd_interp.nii -add T2TRA_thr_bias_corr_reg.nii.gz -div 2 T2TRA_thr_bias_corr_avg.nii.gz" v

			exe "mv T2TRA_thr_bias_corr.nii.gz T2TRA_thr_bias_corr_old.nii.gz" v
			exe "mv T2TRA_thr_bias_corr_avg.nii.gz T2TRA_thr_bias_corr.nii.gz" v

			exe "rm odd_to_even.* slice_reg* T2TRA_thr_bias_corr_odd_interp.* T2TRA_thr_bias_corr_even_interp.* T2TRA_thr_bias_corr_reg.*" v

			show "Slice-by-slice correction of T2TRA image is done." y

		else
			show "Slice-by-slice correction of T2TRA image has been done before." g
		fi
	fi

#----------------------------------------------------------------------------------------------------
# Resampling and cropping of T2 sagittal (T2SAG) (for fitting with resolution of T2TRA image)

	if [ ! -f T2SAG_resample.nii.gz ];then

		show "Starting resampling and cropping of T2SAG image..." y

		# Resample of T2 sagittal image from 0.27 x 0.27 x 1.3 to 0.27 x 0.27 x 0.35

		# Get pixdim1 and pixdim2 of T2SAG
		pixdim1=$(fslval T2SAG_thr_bias_corr.nii.gz pixdim1)
		pixdim2=$(fslval T2SAG_thr_bias_corr.nii.gz pixdim2)
		# Get pixdim1 of T2TRA (will be pixdim3 for T2SAG)
		pixdim3=$(fslval T2TRA_thr_bias_corr.nii.gz pixdim1)


		#exe "sct_resample -i T2SAG_thr_bias_corr.nii.gz -mm 0.279018x0.279018x0.351562 -x spline -o T2SAG_resample.nii.gz" v
		exe "sct_resample -i T2SAG_thr_bias_corr.nii.gz -mm ${pixdim1}x${pixdim2}x${pixdim3} -x spline -o T2SAG_resample.nii.gz" v

		# Thresholding to interval [0 1] (after resampling some negative intensities are included)
		exe "fslmaths T2SAG_resample.nii.gz -thr 0 T2SAG_resample.nii.gz" v

		# Crop T2SAG for faster processing

		while read line; do
			if [[ $line = ${SUB}* ]];then
				START=$(cut -c 6-10 <<< $line);
				END=$(cut -c 12- <<< $line);
			fi;
		done < ${BINDIR}/crop_t2sag.txt

		if [ $START != "" ] && [ $END != "" ];then
			exe "sct_crop_image -i T2SAG_resample.nii.gz -o T2SAG_resample.nii.gz -v 0 -start $START -end $END -dim 1" v
		else
			exe "sct_crop_image -i T2SAG_resample.nii.gz -o T2SAG_resample.nii.gz -v 0 -start 0.25 -end 0.75 -dim 1" v
		fi

		show "Resampling and cropping of T2SAG image is done." y

	else
		show "Resampling and cropping of T2SAG image has been done before." g
	fi

#----------------------------------------------------------------------------------------------------
# Spinal cord segmentation of T2 sagittal image (T2SAG) either by sct_deepseg_sc or sct_propseg function

	if [ ! -f T2SAG_resample_seg.nii.gz ];then

		# SC segmentation using sct_propseg algorithm (older one)
		if [ $DEEPSEG == 0 ];then

			if [ ! -f T2SAG_resample_init_centerline.nii.gz ];then

				show "Starting T2SAG spinal cord segmentation (without using init-centerline)..." y

				# Spinal cord segmentation - modified by JV (added CSF segmentation, change radius to 4 and set max-deformation to default value)
				exe "$SCTDIR/bin/sct_propseg -i T2SAG_resample.nii.gz -c t2 -detect-n 6 -centerline-binary -radius 6 -max-deformation 5 -CSF" v
				#sct_propseg -i T2SAG_resample.nii.gz -c t2 -detect-n 6 -centerline-binary -radius 4 -max-deformation 5 -CSF

				show "T2SAG spinal cord segmentation (without using init-centerline) is done." y

			# Using of manually created init-centerline
			else

				show "Starting T2SAG spinal cord segmentation (using init-centerline)..." y

				# !!!!!!!! Problem has occured with using previous command on patients - it is neccesary to manually create init-centerline, added by JV, 13.11.17
				exe "sct_propseg -i T2SAG_resample.nii.gz -c t2 -detect-n 6 -centerline-binary -radius 4 -max-deformation 5 -CSF -init-centerline T2SAG_resample_init_centerline.nii.gz" v

				show "T2SAG spinal cord segmentation (using init-centerline) is done." y

			fi

		# SC segmentation using sct_deepseg algortim based on CNN (newer one)
		elif [ $DEEPSEG == 1 ];then
			show "Starting T2SAG spinal cord segmentation using sct_deepseg_sc function." y

			exe "sct_maths -i T2SAG_resample.nii.gz -mul 1000 -o T2SAG_resample_scaled.nii.gz" v
			exe "sct_deepseg_sc -i T2SAG_resample_scaled.nii.gz -c t2" v
			exe "mv T2SAG_resample_scaled_seg.nii.gz T2SAG_resample_seg.nii.gz" v
			exe "rm T2SAG_resample_scaled.nii.gz" v

			show "T2SAG spinal cord segmentation using sct_deepseg_sc is done." y
		fi

	else
		show "T2SAG spinal cord segmentation has been done before." g
	fi

#----------------------------------------------------------------------------------------------------
# While loop for opening of T2 sagittal image (T2SAG) and manual enetering of position of C2/3 level

# Control if text file vert_index.txt contains number of positon C2/C3 slice, if not script ask you for enter it. JV, 28.6.2018

# check if vert_index.txt contains subjectID
	# txt file does not containt subID
	if [[ `cat ${BINDIR}/vert_index.txt | grep $SUB` = "" ]];then

		show "Textfile ${BINDIR}/vert_index.txt does not contain postion of C2/C3 slice of current subject!" w
		show "fsleyes ${DATA}/${SUB}/Results/Anat_Preproc/T2SAG_resample.nii.gz -dr 0 0.5 ${DATA}/${SUB}/Results/Anat_Preproc/T2SAG_resample_seg.nii.gz -cm red -a 50 &"
		show "Open FSLeyes by previous command, find C2/C3 slice and type subjectID and y voxel location and press [ENTER] (eg. ${SUB} 456) or if the SC segmentation is bad type bad: "

		# Read from answer from command line
		while read;do
			# bad - segmentation is totaly wrong
			if [[ ${REPLY} = bad ]];then
				echo $SUB $REPLY >> ${BINDIR}/vert_index.txt

				BAD_T2SAG_SEG=1;	# if T2SAG SC segmentation is bad (e.g. due to poor contrast), set this variable to 1 for skip T2SAG labeling, T2SAG registration to template and T2 fusing

				show "$SUB $REPLY" v
				show "T2SAG SC segmentation is bad. Skipping T2SAG vert. labeling and registration to template!!!" w
				break
			# incorrect input - try again
			elif [[ ${#REPLY} -ne 9 ]];then
				case $REPLY in
					* ) show "Incorrect input, enter it again (eg. $SUB 456): " ;;
				esac
			# save answer to file
			else
				echo $REPLY >> ${BINDIR}/vert_index.txt
				show "$REPLY" v
				break
			fi
		done

	# txt file contains subID and T2SAG is bad
	elif [[ `cat ${BINDIR}/vert_index.txt | grep $SUB` =~ bad ]];then
		BAD_T2SAG_SEG=1;
		show "T2SAG SC segmentation is bad. Skipping T2SAG vert. labeling and registration to template!!!" w
	fi

#----------------------------------------------------------------------------------------------------
# Vertebrae labeling of T2SAG

# In vert_index.txt there are slice numbers for every subject
	# get C2/3 slice from txt file
	while read line; do
		if [[ $line = ${SUB}* ]];then
			VERT=$(cut -c 6- <<< $line);
		fi;
	done < ${BINDIR}/vert_index.txt


	# Do only if T2SAG is not wrong (if T2SAG is wrong it is skipped)
	if [[ $BAD_T2SAG_SEG = 0 ]];then
		if [ ! -f T2SAG_resample_seg_labeled.nii.gz ];then

			show "Starting T2SAG vertebrae labeling..." y

			# Vertebrae labeling
			exe "sct_label_vertebrae -i T2SAG_resample.nii.gz -s T2SAG_resample_seg.nii.gz -c t2 -initz ${VERT},2" v
			# Create labels at C2 and C7 vertebral levels
			exe "sct_label_utils -i T2SAG_resample_seg_labeled.nii.gz -vert-body 3,7" v

			show "T2SAG vertebrae labeling is done." y

		else
			show "T2SAG vertebrae labeling has been done before." g
		fi
#----------------------------------------------------------------------------------------------------
# Register T2SAG to template (T2SAG SC segmentation is necessary)

		if [ $SUB != 2648B ];then		# This subject has bad T2SAG FOV
			if [ ! -f warp_template2anat.nii.gz ];then

				show "Starting T2SAG registration to template..." y

				exe "sct_register_to_template -i T2SAG_resample.nii.gz -s T2SAG_resample_seg.nii.gz -l labels.nii.gz -c t2" v
				#fsleyes /usr/local/lib/sct/data/PAM50/template/PAM50_t2.nii.gz anat2template.nii.gz &

				# Warp templates without the white matter atlas (we don't need it at this point)
				#exe "sct_warp_template -d T2SAG_resample.nii.gz -w warp_template2anat.nii.gz -a 0" v	#Not necesary. Template and atlas will be warped to T2TRA below

				show "T2SAG registration to template is done." y

			else
				show "T2SAG registration to template has been done before." g
			fi
		fi
	fi

#----------------------------------------------------------------------------------------------------
# Spinal cord segmentation of T2TRA etiher by sct_propseg or using T2SAGseg as init-centerline or manually created init-centerline) or by sct_deepseg_sc

	if [ ! -f T2TRA_thr_bias_corr_seg.nii.gz ];then

		show "Starting T2TRA spinal cord segmentation..." y

		exe "sct_register_multimodal -i T2SAG_resample_seg.nii.gz -d T2TRA_thr_bias_corr.nii.gz -identity 1 -x nn" v	# necessary for labeling or as init-centerline in sct_propseg

		# T2TRA SC segmentation using sct_propseg algortihm
		if [ $DEEPSEG == 0 ];then

			# Segment T2TRA SC using T2SAG segmentation as init-centerline
			if [ ! -f T2TRA_thr_bias_corr_init_centerline.nii.gz ];then
				# Settings for healthy controls
				if [ ${DATA: -6} == NA-CSD ];then
					exe "sct_propseg -i T2TRA_thr_bias_corr.nii.gz -c t2 -init-centerline T2SAG_resample_seg_reg.nii.gz -detect-n 6 -centerline-binary -radius 6 -max-deformation 4 -CSF" v

					show "T2TRA spinal cord segmentation of healthy subject (using T2SAGseg as initcenterline) is done." y
				# Settings for patients with compression
				elif [ ${DATA: -6} == MA-CSD ];then
					exe "sct_propseg -i T2TRA_thr_bias_corr.nii.gz -c t2 -init-centerline T2SAG_resample_seg_reg.nii.gz -detect-n 6 -centerline-binary -radius 4 -max-deformation 4 -CSF" v

					show "T2TRA spinal cord segmentation of pacient (using T2SAGseg as initcenterline) is done." y

				fi
			# Segment T2TRA SC using manually created centerline
			else
				exe "sct_propseg -i T2TRA_thr_bias_corr.nii.gz -c t2 -init-centerline T2TRA_thr_bias_corr_init_centerline.nii.gz -detect-n 6 -centerline-binary -radius 4 -max-deformation 6 -CSF" v

				show "T2TRA spinal cord segmentation (using manually created init-centerline) is done." y

			fi

		# T2TRA SC segmentation using sct_deepseg algorithm (based on CNN)
		elif [ $DEEPSEG == 1 ];then

			exe "sct_maths -i T2TRA_thr_bias_corr.nii.gz -mul 1000 -o T2TRA_thr_bias_corr_scaled.nii.gz" v
			exe "sct_deepseg_sc -i T2TRA_thr_bias_corr_scaled.nii.gz -c t2" v
			exe "mv T2TRA_thr_bias_corr_scaled_seg.nii.gz T2TRA_thr_bias_corr_seg.nii.gz" v
			exe "rm T2TRA_thr_bias_corr_scaled.nii.gz" v

			# Few subjects with bad results
			if [ $SUB == 2246B ] || [ $SUB == 2284B ];then
				show "Automacic SC segmentation is probably wrong - need to be manually corected using FSLeyes!!!" w
				exit
			fi

			show "T2TRA spinal cord segmentation using sct_deepseg_sc is done." y
		fi

	else
		show "T2TRA spinal cord segmentation has been done before." g
	fi

	#fsleyes T2TRA_thr_bias_corr.nii.gz T2TRA_thr_bias_corr_seg.nii.gz -cm red  &

#----------------------------------------------------------------------------------------------------
# Segment GM from T2TRA either by sct_segment_graymatter or by sct_deepseg_gm
# GM segmentation will be checked using FSLeyes with T2TRA labeling

	if [ ! -f T2TRA_thr_bias_corr_gmseg.nii.gz ];then

		show "Starting GM segmentation of T2TRA image..." y

		# Segment GM and WM using sct_segment_graymatter algorithm (older one)
		if [ $DEEPSEG == 0 ];then

			exe "sct_segment_graymatter -i T2TRA_thr_bias_corr.nii.gz -s T2TRA_thr_bias_corr_seg.nii.gz" v

			show "GM and WM segmentation of T2TRA image is done." y

		# Segment GM using sct_deepseg_gm algortithm (newer on based on CNN)
		elif [ $DEEPSEG == 1 ];then
			exe "sct_maths -i T2TRA_thr_bias_corr.nii.gz -mul 1000 -o T2TRA_thr_bias_corr_scaled.nii.gz" v
			exe "sct_deepseg_gm -i T2TRA_thr_bias_corr_scaled.nii.gz" v
			exe "mv T2TRA_thr_bias_corr_scaled_gmseg.nii.gz T2TRA_thr_bias_corr_gmseg.nii.gz" v
			exe "rm T2TRA_thr_bias_corr_scaled.nii.gz" v

			show "T2TRA GM segmentation using sct_deepseg_gm is done." y
		fi

	else
		show "GM segmentation of T2TRA image has been done before." g
	fi

	#fsleyes T2TRA_thr_bias_corr.nii.gz ./T2TRA_thr_bias_corr_gmseg.nii.gz -cm red-yellow ./T2TRA_thr_bias_corr_wmseg.nii.gz -cm blue-lightblue &

#----------------------------------------------------------------------------------------------------
# Renaming of manually corrected segmentations

	# SC
	if [ -f T2TRA_thr_bias_corr_seg_TH.nii.gz ] && [ ! -f T2TRA_thr_bias_corr_seg_deepseg.nii.gz ];then

		show "Script has found manually corrected SC segmentation. This segmentation will be used in next analysis." w
		exe "mv T2TRA_thr_bias_corr_seg.nii.gz T2TRA_thr_bias_corr_seg_deepseg.nii.gz" v
		exe "mv T2TRA_thr_bias_corr_seg_TH.nii.gz T2TRA_thr_bias_corr_seg.nii.gz" v

		echo "$SUB SC" >> ${BINDIR}/CSD_manually_corrected_segmentation.txt

	fi


	# GM
	if [ -f T2TRA_thr_bias_corr_gmseg_TH.nii.gz ] && [ ! -f T2TRA_thr_bias_corr_gmseg_deepseg.nii.gz ];then

		show "Script has found manually corrected GM segmentation. This segmentation will be used in next analysis." w
		exe "mv T2TRA_thr_bias_corr_gmseg.nii.gz T2TRA_thr_bias_corr_gmseg_deepseg.nii.gz" v
		exe "mv T2TRA_thr_bias_corr_gmseg_TH.nii.gz T2TRA_thr_bias_corr_gmseg.nii.gz" v

		echo "$SUB GM" >> ${BINDIR}/CSD_manually_corrected_segmentation.txt

	fi

#----------------------------------------------------------------------------------------------------
# While loop for spinal cord labeling of T2TRA image
# this part of code calls t2tra_labeling function defined at the beginning of the script

	if [ ! -f T2TRA_thr_bias_corr_seg_labeled.nii.gz ];then

		show "Starting T2TRA vertebrae labeling..." y

		# Call function for labeling
		t2tra_labeling

		# Backup original labeling
		cp T2TRA_thr_bias_corr_seg_labeled.nii.gz T2TRA_thr_bias_corr_seg_labeled_orig.nii.gz

		# Loop for checking if labeling is correct or not. If YES -> continue, If NOT -> repeat it with entered parameters
		while read;do
			# Incorrect input
			if [[ ${#REPLY} -ne 11 ]];then
				case $REPLY in
					[yYaA]* ) break ;;
					* ) show "Incorrect input, enter it again (eg. y/Y/a/A or $SUB 45 10): " ;;
				esac
			# Correct input
			else
				if [[ `cat ${BINDIR}/vert_index_tra.txt | grep $SUB` = "" ]];then	# Check if textfile contains subjectID or not
					echo $REPLY >> ${BINDIR}/vert_index_tra.txt			# Add line into textfile
				else
					sed -i "/${SUB}/d" ${BINDIR}/vert_index_tra.txt			# Delete line if is contained already
					echo $REPLY >> ${BINDIR}/vert_index_tra.txt			# Add line into textfile
					show "$REPLY" v						# Save answer into log
				fi

				t2tra_labeling								# Call function for labeling
			fi
		done

		# WM masking (Subtraction SC and GM). It is in this step because in previous FSLeyes user can correct GM segmentation if it is wrong.
		exe "fslmaths T2TRA_thr_bias_corr_seg.nii.gz -sub T2TRA_thr_bias_corr_gmseg.nii.gz -thr 0 T2TRA_thr_bias_corr_wmseg.nii.gz" v

		show "T2TRA vertebrae labeling is done." y

	else
		show "T2TRA vertebrae labeling has been done before." g
	fi

#----------------------------------------------------------------------------------------------------
# Register template and atlas (PAM50) to T2 transversal image (T2TRA space)

	if [ ! -f warp_template2T2TRA.nii.gz ];then

		show "Starting template and atlas registration to T2TRA..." y

		# Register template to T2TRA
		# T2SAG SC segmentation is ok and warping field from prevous template to T2SAG registration is now used as initwarp
		if [[ $BAD_T2SAG_SEG = 0 ]];then
			exe "sct_register_multimodal -i $SCTDIR/data/PAM50/template/PAM50_t2s.nii.gz -iseg $SCTDIR/data/PAM50/template/PAM50_cord.nii.gz -d T2TRA_thr_bias_corr.nii.gz -dseg T2TRA_thr_bias_corr_seg.nii.gz -param step=1,type=seg,algo=centermassrot,smooth=2:step=2,type=seg,algo=bsplinesyn,slicewise=1,iter=3,smooth=1 -initwarp warp_template2anat.nii.gz -owarp warp_template2T2TRA.nii.gz" v
		# T2SAG SC is bad and registration do NOT use any initwarp
		else
			show "Script detected bad T2SAG SC segmentation. Template and atlas registration will be performed without init warpfield!!!" w
			exe "sct_register_multimodal -i $SCTDIR/data/PAM50/template/PAM50_t2s.nii.gz -iseg $SCTDIR/data/PAM50/template/PAM50_cord.nii.gz -d T2TRA_thr_bias_corr.nii.gz -dseg T2TRA_thr_bias_corr_seg.nii.gz -param step=1,type=seg,algo=centermassrot,smooth=2:step=2,type=seg,algo=bsplinesyn,slicewise=1,iter=3,smooth=1 -owarp warp_template2T2TRA.nii.gz" v
		fi

		# Warp template and atlas (PAM50) to T2TRA - warped folder label/ will be used in sct_register_graymatter command
		exe "sct_warp_template -d T2TRA_thr_bias_corr.nii.gz -w warp_template2T2TRA.nii.gz -a 0" v

		show "Template and atlas registration to T2TRA is done." y

	else
		show "Template and atlas registration to T2TRA has been done before." g
	fi
	#fsleyes T2TRA_thr_bias_corr.nii.gz PAM50_t2s_reg.nii.gz &

#----------------------------------------------------------------------------------------------------
# Improve template registration by acounting WM/GM shape (requires correct WM and GM segmentation in T2TRA space)

	if [ ! -f warp_template2T2TRA_reg_gm.nii.gz ];then

		show "Starting template GM and WM registration to GM and WM of T2TRA image..." y

		# Improve template registration by accounting GM and WM shape
		exe "sct_register_graymatter -gm T2TRA_thr_bias_corr_gmseg.nii.gz -wm T2TRA_thr_bias_corr_wmseg.nii.gz -w warp_template2T2TRA.nii.gz" v
		# Warp template and atlas (PAM50) (this time corrected for internal structure)
		exe "sct_warp_template -d T2TRA_thr_bias_corr.nii.gz -w warp_template2T2TRA_reg_gm.nii.gz -a 1" v

		show "Improving of template and atlas registration by accounting GM and WM shape is done." y

	else
		show "Improving of template and atlas registration by accounting GM and WM shape has been done before." g
	fi

# fsleyes PAM50_t2s_reg.nii.gz T2TRA_thr_bias_corr.nii.gz label/template/PAM50_cord.nii.gz label/template/PAM50_levels.nii.gz

}
# end of anat_analysis function

#------------------------

# Call script for diffusion analysis
diff_analysis()
{
	/md1/ing/bin/NeuroImage/diff_analysis.sh $DATA $SUB $SEQORDER
}

#------------------------

init_task $@
