#!/bin/bash

##
## for S in `ls -d ???_?_20??????`; do echo " ++++ $S/DICOMS/"; find $S/DICOMS -maxdepth 3 -type d -name "*T1w*" -print ; done
#  ++++ 001_A_20160601/DICOMS/
# 001_A_20160601/DICOMS/PUT1w_BRAVO_sag_3
#  ++++ 001_C_20170531/DICOMS/
# 001_C_20170531/DICOMS/PUT1w_BRAVO_sag_300
# 001_C_20170531/DICOMS/T1w_BRAVO_sag_3
#  ++++ 002_A_20160608/DICOMS/
# 002_A_20160608/DICOMS/PUT1w_BRAVO_sag_3
#  ++++ 003_A_20160608/DICOMS/
# 003_A_20160608/DICOMS/PUT1w_BRAVO_sag_3
#  ++++ 003_C_20200219/DICOMS/
# 003_C_20200219/DICOMS/ORIG_T1w_BRAVO_sag_3
#  ++++ 003_D_20200311/DICOMS/
# 003_D_20200311/DICOMS/ORIG_T1w_BRAVO_sag_3
#  ++++ 004_A_20160608/DICOMS/
# 004_A_20160608/DICOMS/PUT1w_BRAVO_sag_3
#  ++++ 004_C_20170614/DICOMS/
# 004_C_20170614/DICOMS/PUT1w_BRAVO_sag_600
# 004_C_20170614/DICOMS/T1w_BRAVO_sag_6
#  ++++ 005_A_20160615/DICOMS/
# 005_A_20160615/DICOMS/PUT1w_BRAVO_sag_5
#  ++++ 005_C_20170607/DICOMS/
# 005_C_20170607/DICOMS/PUT1w_BRAVO_sag_600
# 005_C_20170607/DICOMS/T1w_BRAVO_sag_6
#  ++++ 005_E_20190426/DICOMS/
# 005_E_20190426/DICOMS/T1w_BRAVO_sag_4
# 005_E_20190426/DICOMS/ORIG_T1w_BRAVO_sag_40004



Usage() {
	echo
	echo "Usage: `basename $0` <ssid>"
	echo
	echo "Example: `basename $0` NSEPP001"
	echo
	exit 1
}


SCRIPT=$(python -c "from os.path import abspath; print(abspath('$0'))")
SCRIPTSDIR=$(dirname $SCRIPT)
if [[ -z "${RAW_DATA_DIR}" ]]; then
	source "$SCRIPTSDIR/FORBOW_SetupEnvironment.sh"
fi
FORCE_OVERWRITE="no"
if [ "$1" == "-f" ]; then
	FORCE_OVERWRITE="yes"
	shift ;
fi

[ "$#" -lt 1 ] && Usage

## SPECIFIC VARIABLES
DO_DEFACING="yes"
DO_CLEANUP="no"
DEFACE_METHOD="ants"		##options=[ants,fsl]
DO_BASIC_QC="yes"
MNI_BN="$MNI152DIR/MNI152"
MNI_T1_1mm="${MNI_BN}_T1_1mm.nii.gz"
MNI_T1_1mm="${MNI_BN}_T1_1mm.nii.gz"
MNI_T1="${MNI_BN}_T1_1mm.nii.gz"
MNI_T2="${MNI_BN}_T2_1mm.nii.gz"
FACEMASK_Small="${MNI_BN}_T1_1mm_facemask.nii.gz"
FACEMASK_Big="${MNI_BN}_T1_1mm_BigFoV_facemask.nii.gz"
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=4
DWI_peAP_NUMVOLS="33"
DWI_pePA_NUMVOLS="8"
REDO_QC="yes"

SUBJECTS="$@"
for Subj in $SUBJECTS ; do 
	
	S=$(basename $Subj)
	if [ "${#S}" -eq 14 ]; then
		SSID="${S:0:5}"	
		SDIR="${RAW_DATA_DIR}/$S"
		S=$SSID
	elif [ "${#S}" -eq 5 ]; then
		SDIR=$(find ${RAW_DATA_DIR} -maxdepth 1 -type d -name "${S}_20??????" -print)
	else
		echo " * ERROR: cannot determine subject from input=[$Subj]"
		continue
	fi
	SDIR_BN=$(basename $SDIR)
	if [ "${#SDIR_BN}" -eq 14 ]; then
		SubjAcqDate=${SDIR_BN:6:8}
	else
		echo "*** ERROR: failed to determine AcqDate from SDIR=${SDIR}..."
		continue
	fi
	if [ ! -d "$SDIR" ]; then
		echo "*** ERROR: cannot find rawdata directory for subject=[$S] in ${RAW_DATA_DIR}/"
		continue
	fi
	niiDIR="$SDIR/NIFTIS"
	cd $niiDIR/
	echo "--------- running `basename $0` for `pwd`/, starting on: `date`"
	
	if [ "$DO_DEFACING" == "yes" ]; then
		echo " ++ Defacing Anats"
		mkdir -p $niiDIR/1_defacing_work/ && cd $niiDIR/1_defacing_work/
		for seq in t1w_bravo_orig  t2prep_flair_promo_orig ; do 
			img="${S}_${seq}"
			if [ ! -r "../${img}.nii.gz" ]; then 
				continue
			fi
			ln -sf ../${img}.nii.gz ${seq}.nii.gz
			if [ ! -r "${seq}_n4.nii.gz" ]; then
				## Do quick Bias Correction
				${ANTSPATH}/N4BiasFieldCorrection -d 3 -i ${seq}.nii.gz -o ${seq}_n4.nii.gz -s 4 -b [200] -c [50x50x50x50,0.000001]
			fi
			MNI="${MNI_BN}_T1"
			isT2=$(echo "$seq" | grep -i "t2")
			if [ "$seq" == "$isT2" ]; then
				MNI="${MNI_BN}_T2"
			fi
			if [ ! -r "../${img}_defaced.nii.gz" ]; then
				if [ "$DEFACE_METHOD" == "ants" ]; then				
					ST=$SECONDS;
					## Do quick registration to MNI152 with ANTs
					${ANTSPATH}/antsRegistrationSyNQuick.sh -d 3 -f ${MNI}_1mm.nii.gz -m ${seq}_n4.nii.gz -o ${seq}_n4_2mni_ -t s
					## Apply registration to move mask from MNI -> img
					${ANTSPATH}/antsApplyTransforms \
						--dimensionality 3 \
						--interpolation GenericLabel \
						--reference-image ${seq}_n4.nii.gz \
						--input $FACEMASK_Big \
						--output ../${img}_facemask.nii.gz \
						--transform ${seq}_n4_2mni_1InverseWarp.nii.gz \
						--transform [ ${seq}_n4_2mni_0GenericAffine.mat , 1 ] \
						--verbose 0
					fslmaths ../${img}_facemask.nii.gz -binv -mul ${seq}.nii.gz ../${img}_defaced.nii.gz
					echo " + finished ants defacing on ${seq}.nii.gz -> ${MNI} in $(($SECONDS-$ST)) secs"
				elif [ "$DEFACE_METHOD" == "fsl" ]; then
					ST=$SECONDS;
					echo " ++ running initial flirt to register ${seq}_n4.nii --> ${MNI}_2mm.nii"
					flirt -in ${seq}_n4.nii -ref ${MNI}_2mm.nii -dof 12 -coarsesearch 30 -omat ${seq}_2mni2mm.mat -o ${seq}_2mni2mm.nii
					echo " ++ running final flirt to register ${MNI}_1mm.nii --> ${seq}_n4.nii, using -init ${seq}_2mni2mm.mat"
					flirt -in ${seq}_n4.nii -ref ${MNI}_1mm.nii -dof 12 -init ${seq}_2mni2mm.mat -omat ${seq}_2mni1mm.mat -o ${seq}_2mni1mm.nii
					convert_xfm -omat ${seq}_2mni1mm_BigFoV.mat -concat $FSLDIR/etc/flirtsch/MNI_to_MNI_BigFoV_facemask.mat ${seq}_2mni1mm.mat
					convert_xfm -omat ${seq}_2mni1mm_BigFoV_inverse.mat -inverse  ${seq}_2mni1mm_BigFoV.mat
					echo " ++ registering $FACEMASK_Big --> ${seq}"
					flirt -in $FACEMASK_Big -ref ${seq} -init ${seq}_2mni1mm_BigFoV_inverse.mat -applyxfm -interp nearestneighbour -o ../${img}_facemask.nii
					echo " ++ masking ${seq} with ../${img}_facemask.nii.gz"
					fslmaths ../${img}_facemask.nii.gz -binv -mul ${seq}.nii.gz ../${img}_defaced.nii.gz
					echo " + finished fsl defacing ${seq}.nii.gz, in $(($SECONDS-$ST)) secs"
				else
					echo "*** ERROR: variable 'DEFACE_METHOD' must be set = [ants,fsl]."
					exit 3
				fi
			fi	
			if [ "$DO_CLEANUP" == "yes" ]; then
				cd $niiDIR/
				rm -rf $niiDIR/1_defacing_work/
			fi
		done
		echo " ++ finished Defacing Anats"
	fi
	
	if [ "$DO_BASIC_QC" == "yes" ]; then
		echo " --- Creating Basic QC reports"
		qcDIR="$niiDIR/QC"
		if [ "$REDO_QC" == "yes" ]; then
			rm -rf ${qcDIR}/ >/dev/null
		fi
		mkdir -p $qcDIR/
		qcHTML="$qcDIR/_QC_report.html"
		mkdir -p $niiDIR/2_qc_work/ && cd $niiDIR/2_qc_work/ 
		echo "<HTML><TITLE></TITLE><BODY BGCOLOR=\"#ffffff\"><h1>QC rawdata report for ${S}<h1>" >$qcHTML
		echo "<br><hr><h2 style=”color: red; font-weight: bold;”>$SDIR/</h2>" >>$qcHTML
		##for seq in t1w_bravo_orig t1w_bravo_pure t2prep_flair_promo_orig t2w_cube_promo_orig t2w_cube_promo_pure t2w_cube_orig t2w_cube_pure dwi_peAP dwi_pePA rs_SBRef rs_peAP rs_pePA ; do
		for seq in t1w_bravo_orig t2prep_flair_promo_orig dwi_peAP dwi_pePA rs_SBRef rs_peAP rs_pePA ; do
			img="${S}_${seq}"
			ln -sf ../${img}.nii.gz ./${seq}.nii.gz
			if [ "${seq:0:3}" == "dwi" ]; then
				if [ ! -r "${seq}_nodif.nii.gz" ]; then
					echo " +-+ preparing seq=$img for QC..."
					bvals="../${img}.bval"
					B0Vols=$(python -c "x=[str(i) for i,v in enumerate(open('$bvals','r').read().split()) if float(v) < 50.0]; print(','.join(x))")
					fslselectvols -i ${seq} -o ${seq}_B0s --vols=${B0Vols}
					mcflirt -in ${seq}_B0s -meanvol -spline_final -o ${seq}_B0s_mc
					fslmaths ${seq}_B0s_mc -thr 0 -Tmean ${seq}_nodif
				fi
				img="${seq}_nodif"
			elif [ "${seq:0:2}" == "rs" ]; then 
				if [ ! -r "${seq}_mean.nii.gz" ]; then
					echo " +-+ preparing seq=$img for QC..."
					mcflirt -in ${seq}.nii.gz -meanvol -spline_final -o ${seq}_mc
					fslmaths ${seq}_mc -thr 0 -Tmean ${seq}_mean
				fi
				img="${seq}_mean"
			else
				img="${seq}"
			fi
			if [ -f "${seq}.nii.gz" -a ! -f "$qcDIR/${seq}_9x1.png" ]; then
				#echo " ++ doing QC on $img"
				create_9x1_slicer ${img}.nii.gz $qcDIR/${seq}_9x1.png
				hdr=$(fslinfo $niiDIR/${S}_${seq}.nii.gz | tr -s '\n' '\t' | awk '{printf("matrix=%dx%dx%d@%s, nvols=%d, pixdims=%.2fx%.2fx%.2f,TR=%.3f\n",$4,$6,$8,$2,$10,$14,$16,$18,$20)}')
				echo "<p><strong>${seq}  ($hdr):</strong><br><img src=\"${seq}_9x1.png\" WIDTH=2000></p>" >>$qcHTML
			else
				cp -p ${SCRIPTSDIR}/missing_overlay_9x1.png $qcDIR/${seq}_9x1.png
				echo "<p><strong>${seq}  ():</strong><br><img src=\"${seq}_9x1.png\" WIDTH=2000></p>" >>$qcHTML
				echo -e "\n *** WARNING: missing sequence = $niiDIR/${img}.nii.gz"
			fi
			if [ -r "${img}.nii.gz" -a "${seq:0:3}" != "dwi" -a "${seq:0:2}" != "rs" ]; then
				echo " ++ creating a quick overlay to visualize facemaks on ${seq}"
				t1w98P=$(fslstats ${img}.nii.gz -P 98 | LC_ALL=C xargs /usr/bin/printf "%.0f\n")
				$FSLDIR/bin/overlay 1 0 ${img}.nii.gz 0 $t1w98P $niiDIR/${S}_${seq}_facemask.nii.gz 1 2 ${img}_facemask_rendered.nii.gz
				create_8x3_overlay ${img}_facemask_rendered.nii.gz $qcDIR/${seq}_facemask_overlay_8x3.png ${FSLDIR}/etc/luts/render1t.lut
				echo "<p><strong>${seq}_facemask.nii </strong><br><img src=\"${seq}_facemask_overlay_8x3.png\" WIDTH=2000></p>" >>$qcHTML
			fi
		done
		echo "<br></BODY></HTML>" >>$qcHTML
		if [ "$DO_CLEANUP" == "yes" ]; then
			cd $niiDIR/
			rm -rf $niiDIR/2_qc_work/
		fi
	fi
	
	echo "--------- finished `basename $0` on `pwd`, at: `date`"	
	ls -l 
done

