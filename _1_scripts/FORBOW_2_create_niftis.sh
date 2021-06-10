#!/bin/bash
#
# SCRIPTNAME = rawdata/_1_scripts/FORBOW_2_create_niftis_20201116.sh
# 
# This script does three things: 
#  1) converts T1w/T2w/DWI /DICOMS/ into /NIFTIS/ using dcm2niix_v20201102. Multiple T1w and 
#     T2w sequences where naming on console changed several times, each is converted and named 
#     as follows:
#
#    T1w:
#              ** before 2017/12/01 **
#      /DICOMS/T1w_BRAVO_sag/		-> /NIFTIS/ssid_t1w_bravo_orig.nii.gz
#      /DICOMS/PUT1w_BRAVO_sag/		-> /NIFTIS/ssid_t1w_bravo_pure.nii.gz
#              ** after 2017/12/01 **
#      /DICOMS/ORIG_T1w_BRAVO_sag/	-> /NIFTIS/ssid_t1w_bravo_orig.nii.gz
#      /DICOMS/T1w_BRAVO_sag/ 		-> /NIFTIS/ssid_t1w_bravo_pure.nii.gz
#
#    T2w-CUBE:
#              ** before 2017/12/01 **
#      /DICOMS/T2w_CUBE_sag/		-> /NIFTIS/ssid_t2w_cube_orig.nii.gz
#      /DICOMS/PUT2w_CUBE_sag/ 		-> /NIFTIS/ssid_t2w_cube_pure.nii.gz
#              ** after 2017/12/01 **
#      /DICOMS/ORIG_T2w_CUBE_sag/ 	-> /NIFTIS/ssid_t2w_cube_orig.nii.gz
#      /DICOMS/T2w_CUBE_sag/ 		-> /NIFTIS/ssid_t2w_cube_pure.nii.gz
#
#    T2w-CUBE-PROMO:
#              ** before 2017/12/01 **
#      /DICOMS/T2w_CUBE_PROMO_sag/		-> /NIFTIS/ssid_t2w_cube_orig.nii.gz
#              ** after 2017/12/01 **
#      /DICOMS/ORIG_T2w_CUBE_PROMO_sag/ -> /NIFTIS/ssid_t2w_cube_promo_orig.nii.gz
#
#    T2wFLAIR_Prep-PROMO:
#      /DICOMS/T2PREP_FLAIR_PROMO/	-> /NIFTIS/ssid_t2prep_flair_promo_orig.nii.gz
#
#    DWI:
#      /DICOMS/DTI_30_DIR+3B0/		-> /NIFTIS/ssid_dwi_peAP.nii.gz
#      /DICOMS/DTI_Bomap__Rev/		-> /NIFTIS/ssid_dwi_pePA.nii.gz
#
#
#  2) Specifically for /NIHPD_Analysis/, the T1w and T2w sequences below are copied as:
#     /NIFTIS/ssid_t1w_bravo_orig.nii.gz			-> /NIFTIS/ssid_t1.nii.gz
#     /NIFTIS/ssid_t2prep_flair_promo_orig.nii.gz	-> /NIFTIS/ssid_t2.nii.gz
#
#  3) Three RS-Pfile NIFTIS are copied from the /RS/ folder into /NIFTIS/ as follows:
#      /RS/*.nii.gz with vol=500, TR=.950	-> /NIFTIS/ssid_rs_peAP.nii.gz
#      /RS/*.nii.gz with vol=4, TR=.950		-> /NIFTIS/ssid_rs_pePA.nii.gz
#      /RS/*.nii.gz with vol=6, TR=2.80		-> /NIFTIS/ssid_rs_SBRef.nii.gz
#
##########################################################################################


Usage() {
	echo
	echo "Usage: `basename $0` <ssid>"
	echo
	echo "Example: `basename $0` 005_C"
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
dcm2niix=$(which dcm2niix_v20201102) ##use specific version
CONVERT_DICOMS_2_NIFTIS="yes"
ORGANIZE_RS_NIFTIS="yes"
DO_DEFACING="yes"
DO_CLEANUP="yes"
DO_BASIC_QC="yes"

MNI_T1="$MNI152DIR/MNI152_T1_2mm.nii"
MNI_T2="$MNI152DIR/MNI152_T2_2mm.nii"
FACEMASK="$MNI152DIR/MNI152_T1_1mm_facemask.nii"
DWI_peAP_NUMVOLS="8"
DWI_pePA_NUMVOLS="33"
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
	if [ "$FORCE_OVERWRITE" == "yes" ]; then
		echo " * FORCE_OVERWRITE enabled, deleting $niiDIR/"
		rm -rf $niiDIR/ >/dev/null
	fi
	mkdir -p $niiDIR/
	cd $niiDIR/
	echo "--------- running `basename $0` for `pwd`/, starting on: `date`"
	
	REQUIRE_REORIENT_CHECK="no"
	
	## 1) Find and covert T1w & T2w DICOM sequences into NIFTIS
	if [ "$CONVERT_DICOMS_2_NIFTIS" == "yes" ] ; then
		dcmDIR=$(find $SDIR -type d -name "DICOMS" -print | head -1)
		if [ x"$dcmDIR" == "x" ]; then
			echo "*** ERROR: could not locate $SDIR/DICOMS/"
			continue
		fi
		
		### T1w ##########################################################################
		## dicom_hinfo -tag 0020,0011 -tag 0008,0008 005_C_20170607/DICOMS/T1w_BRAVO_sag_6/IM-0003-0001.dcm 		
		# nonPU-orig *	005_C_20170607/DICOMS/T1w_BRAVO_sag_6/IM-0003-0001.dcm 6 ORIGINAL\PRIMARY\OTHER
		# PU-derived -	005_C_20170607/DICOMS/PUT1w_BRAVO_sag_600/IM-0007-0001.dcm 600 DERIVED\PRIMARY\OTHER
		# PU-derived -	005_E_20190426/DICOMS/T1w_BRAVO_sag_4/IM-0001-0001.dcm 4 ORIGINAL\PRIMARY\OTHER
		# nonPU-orig *	005_E_20190426/DICOMS/ORIG_T1w_BRAVO_sag_40004/IM-0005-0001.dcm 40004 ORIGINAL\PRIMARY\OTHER
		##################################################################################
		if [ ! -r "${S}_t1w_bravo_orig.nii.gz" -o ! -r "${S}_t1w_bravo_pure.nii.gz" ]; then
			
			for t1dir in $(find $dcmDIR -type d -name "*T1w_BRAVO_sag*" -print); do
				im1=$(find $t1dir -type f -iname "*.dcm" -print | head -1)
				if [ ! -r "$im1" ]; then
					echo " *** ERROR: cannot find IM-00xx-0001.dcm for T1w sequence=[$t1dir]"
				fi
				bnDIR=$(basename $t1dir)
				isORIG=$(dicom_hinfo -tag 0008,0008 ${im1} | awk '{print $2}')
				serNum=$(dicom_hinfo -tag 0020,0011 ${im1} | awk '{print $2}' | bc)
				acqDate=$(dicom_hinfo -tag 0008,0022 ${im1} | awk '{print $2}')
				#echo " +++ im1=[$im1],[${bnDIR:0:4}] with isORIG=[$isORIG], serNum=[$serNum], acqDate=[$acqDate]"
				if [ "${bnDIR:0:4}" == "ORIG" ]; then
 					echo " + found ORIG T1w sequence=[$t1dir] with isORIG=[$isORIG], serNum=[$serNum], acqDate=[$acqDate]"
 					if [ ! -r "${S}_t1w_bravo_orig.nii.gz" ]; then
						${dcm2niix} -b y -o . -x n -z y -f "${S}_t1w_bravo_orig" $t1dir/
						REQUIRE_REORIENT_CHECK="yes"
					else
						echo "* Warning: found existing file: `pwd`/${S}_t1w_bravo_orig.nii.gz"
					fi
 				elif [ "${bnDIR:0:2}" == "PU" ]; then
 					echo " + found PURE T1w sequence=[$t1dir] with isORIG=[$isORIG], serNum=[$serNum], acqDate=[$acqDate]"
 					if [ ! -r "${S}_t1w_bravo_pure.nii.gz" ]; then
						${dcm2niix} -b y -o . -x n -z y -f "${S}_t1w_bravo_pure" $t1dir/
						REQUIRE_REORIENT_CHECK="yes"
					else
						echo "* Warning: found existing file: `pwd`/${S}_t1w_bravo_pure.nii.gz"
					fi
 				elif [ "${bnDIR:0:3}" == "T1w" -a "${acqDate}" -lt 20171201 ]; then
 					echo " + found ORIG T1w sequence=[$t1dir] with isORIG=[$isORIG], serNum=[$serNum], acqDate=[$acqDate]"
 					if [ ! -r "${S}_t1w_bravo_orig.nii.gz" ]; then
						${dcm2niix} -b y -o . -x n -z y -f "${S}_t1w_bravo_orig" $t1dir/
						REQUIRE_REORIENT_CHECK="yes"
					else
						echo "* Warning: found existing file: `pwd`/${S}_t1w_bravo_orig.nii.gz"
					fi
 				elif [ "${bnDIR:0:3}" == "T1w" -a "${acqDate}" -ge 20171201 ]; then
 					echo " + found PURE T1w sequence=[$t1dir] with isORIG=[$isORIG], serNum=[$serNum], acqDate=[$acqDate]"
 					if [ ! -r "${S}_t1w_bravo_pure.nii.gz" ]; then
						${dcm2niix} -b y -o . -x n -z y -f "${S}_t1w_bravo_pure" $t1dir/
						REQUIRE_REORIENT_CHECK="yes"
					else
						echo "* Warning: found existing file: `pwd`/${S}_t1w_bravo_pure.nii.gz"
					fi
 				else
 					echo " *** unable to determine whether T1w sequence=${t1dir} is ORIGINAL or PURE with isORIG=[$isORIG], serNum=[$serNum], acqDate=[$acqDate]"
 				fi
			done
		fi
		
		### T2w ##########################################################################
		## dicom_hinfo -tag 0020,0011 -tag 0008,0008 005_C_20170607/DICOMS/PUT2w_CUBE_sag_400/IM-0006-0001.dcm		
		# 005_C_20170607/DICOMS/PUT2w_CUBE_sag_400/IM-0006-0001.dcm 400 DERIVED\PRIMARY\OTHER
		# 005_C_20170607/DICOMS/T2PREP_FLAIR_PROMO_5/IM-0002-0001.dcm 5 ORIGINAL\PRIMARY\OTHER
		# 005_C_20170607/DICOMS/T2w_CUBE_sag_4/IM-0001-0001.dcm 4 ORIGINAL\PRIMARY\OTHER
		# 131_D_20190227/DICOMS/ORIG_T2w_CUBE_PROMO_sag_40005/IM-0007-0001.dcm 40005 ORIGINAL\PRIMARY\OTHER
		# 131_D_20190227/DICOMS/T2PREP_FLAIR_PROMO_2/IM-0001-0001.dcm 2 ORIGINAL\PRIMARY\OTHER
		# 131_D_20190227/DICOMS/T2w_CUBE_PROMO_sag_5/IM-0003-0001.dcm 5 ORIGINAL\PRIMARY\OTHER
		##################################################################################
		
		##-------------------------------------------------------------------------------
		## ${S}_t2prep_flair_promo_orig.nii.gz
		if [ ! -r "${S}_t2prep_flair_promo_orig.nii.gz" ]; then
			t2dir=$(find $dcmDIR -type d -name "T2PREP_FLAIR_PROMO*" -print)
			if [ -d "$t2dir" ]; then
				${dcm2niix} -b y -o . -x n -z y -f "${S}_t2prep_flair_promo_orig" $t2dir/
				REQUIRE_REORIENT_CHECK="yes"
			else
				echo "*** ERROR: could not find DICOM sequence 'T2PREP_FLAIR_PROMO'..."
			fi
		fi
		##-------------------------------------------------------------------------------
		## ${S}_t2w_cube_promo_orig.nii.gz
		if [ ! -r "${S}_t2w_cube_promo_orig.nii.gz" ]; then
			## if scanned before 20171201, then ORIG = T2w_CUBE_PROMO
			if [ "$SubjAcqDate" -lt 20171201 ]; then
				t2dir=$(find $dcmDIR -type d -name "T2w_CUBE_PROMO*" -print)
			else
				t2dir=$(find $dcmDIR -type d -name "ORIG_T2w_CUBE_PROMO*" -print)
			fi
			if [ -d "$t2dir" ]; then
				${dcm2niix} -b y -o . -x n -z y -f "${S}_t2w_cube_promo_orig" $t2dir/
				REQUIRE_REORIENT_CHECK="yes"
			else
				echo "* WARNING: could not find DICOM ORIG sequence 'T2w_CUBE_PROMO'"
			fi
		fi
		##-------------------------------------------------------------------------------
		## ${S}_t2w_cube_promo_pure.nii.gz
		if [ ! -r "${S}_t2w_cube_promo_pure.nii.gz" ]; then
			## if scanned before 20171201, then PURE = T2w_CUBE_PROMO
			if [ "$SubjAcqDate" -ge 20171201 ]; then
				t2dir=$(find $dcmDIR -type d -name "T2w_CUBE_PROMO*" -print)
			fi
			if [ -d "$t2dir" ]; then
				${dcm2niix} -b y -o . -x n -z y -f "${S}_t2w_cube_promo_pure" $t2dir/
				REQUIRE_REORIENT_CHECK="yes"
			else
				echo "* Warning: could not find DICOM PURE sequence 'T2w_CUBE_PROMO', collected after 2017/12/01 console update..."
			fi
		fi
		##-------------------------------------------------------------------------------		
		## ${S}_t2w_cube_orig.nii.gz 
		if [ ! -r "${S}_t2w_cube_orig.nii.gz" ]; then
			## if scanned before 20171201, then ORIG = T2w_CUBE
			if [ "$SubjAcqDate" -lt 20171201 ]; then
				t2dir=$(find $dcmDIR -type d -name "T2w_CUBE_sag*" -print)
			else
				t2dir=$(find $dcmDIR -type d -name "ORIG_T2w_CUBE_sag*" -print)
			fi
			if [ -d "$t2dir" ]; then
				${dcm2niix} -b y -o . -x n -z y -f "${S}_t2w_cube_orig" $t2dir/
				REQUIRE_REORIENT_CHECK="yes"
			else
				echo "* Warning: could not find DICOM orig sequence 'T2w_CUBE_sag'..."
			fi
		fi
		##-------------------------------------------------------------------------------		
		## ${S}_t2w_cube_pure.nii.gz 
		if [ ! -r "${S}_t2w_cube_pure.nii.gz" ]; then
			## if scanned before 20171201, then PURE = PUT2w_CUBE_sag
			if [ "$SubjAcqDate" -lt 20171201 ]; then
				t2dir=$(find $dcmDIR -type d -name "PUT2w_CUBE_sag*" -print)
			else
				t2dir=$(find $dcmDIR -type d -name "T2w_CUBE_sag*" -print)
			fi
			if [ -d "$t2dir" ]; then
				${dcm2niix} -b y -o . -x n -z y -f "${S}_t2w_cube_pure" $t2dir/
				REQUIRE_REORIENT_CHECK="yes"
			else
				echo "* Warning: could not find DICOM PURE sequence 'T2w_CUBE_sag'..."
			fi
		fi
		
		
		#### DWI #################################################
		## Following naming scheme to match PE-distortion direction: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup/Faq
		dwi_main="${S}_dwi_pePA"
		if [ ! -r "${dwi_main}.nii.gz" ]; then
			SeqDIR=$(find $dcmDIR -type d -name "DTI_30_DIR+3B0*" -print | tail -1)
			if [ ! -d "$SeqDIR" ]; then
				echo "*** ERROR: cannot locate subject DICOM 'DTI_30_DIR+3B0' directory..."
			else
				echo " + converting $SeqDIR/"
				${dcm2niix} -b y -o . -x n -z y -f "${dwi_main}" $SeqDIR/
				bvalNV=$(cat -v ${dwi_main}.bval | wc -w | bc)
				## verify expected number of volumes or exit with error
				if [ "`fslnvols ${dwi_main}.nii.gz | bc`" -ne "$DWI_pePA_NUMVOLS" -o "${bvalNV}" -ne "$DWI_pePA_NUMVOLS" ]; then
					echo "*** ERROR: ${dwi_main}.nii.gz or ${dwi_main}.bval does not have expected numVols=$DWI_pePA_NUMVOLS"
					exit 3
				fi
				REQUIRE_REORIENT_CHECK="yes"
			fi
		fi
		dwi_rev="${S}_dwi_peAP"
		if [ ! -r "${dwi_rev}.nii.gz" ]; then
			SeqDIR=$(find $dcmDIR -type d -name "DTI_Bomap__Rev*" -print | tail -1)
			if [ ! -d "$SeqDIR" ]; then
				echo "*** ERROR: cannot locate subject DICOM 'DTI_Bomap__Rev' directory..."
			else
				echo " + converting $SeqDIR/"
				${dcm2niix} -b y -o . -x n -z y -f "${dwi_rev}" $SeqDIR/
				bvalNV=$(cat -v ${dwi_rev}.bval | wc -w | bc)
				## verify expected number of volumes or exit with error
				if [ "`fslnvols ${dwi_rev}.nii.gz | bc`" -ne "$DWI_peAP_NUMVOLS" -o "${bvalNV}" -ne "$DWI_peAP_NUMVOLS" ]; then
					echo "*** ERROR: ${dwi_rev}.nii.gz or ${dwi_rev}.bval does not have expected numVols=$DWI_peAP_NUMVOLS"
					exit 3
				fi
				REQUIRE_REORIENT_CHECK="yes"
			fi
		fi
	fi	
	
	
	if [ "$ORGANIZE_RS_NIFTIS" == "yes" ] ; then
		### RS-fMRI ################################################################
		pfileDIR=$(find $SDIR -maxdepth 1 -type d -name "RS" -print | head -1)
		if [ x"$pfileDIR" == "x" ]; then
			echo "*** ERROR: could not locate $SDIR/RS/"
		else
			# -- SBref -->
			if [ ! -r "${S}_rs_SBRef.nii.gz" ]; then
				for p in $(imglob $pfileDIR/*.nii.gz); do
					if [ -r "${p}.nii" ]; then
						pigz -v -b 512 -n -f -6 ${p}.nii
					fi
					H=$(3dinfo -history ${p}.nii.gz | grep '3dLRflip')
					if [ x"$H" == 'x' -a "$SubjAcqDate" -ge 20180101 ]; then
						echo "*** ERROR: RS-NII file=[$p] does NOT have required history containing '3dLRflip'"
						exit 2					
					elif [ x"$H" != 'x' -a "$SubjAcqDate" -lt 20180101 ]; then
						echo "*** ERROR: RS-NII file=[$p] incorrectly has a history containing '3dLRflip'"
						exit 2
					fi
					hdr=$(fslhd $p)
					NV=$(echo "$hdr" | grep '^dim4' | awk {'print $2'} | bc)
					TR=$(echo "$hdr" | grep '^pixdim4' | awk {'print $2'} | bc)
					#echo "-- found P=$p, tr=[$TR],nv=[$NV]"
					if [ "$NV" -eq 6 -a "$TR" == "2.800000" ]; then
						echo "-- copying SBRef RS-fMRI run = $p"
						$FSLDIR/bin/imcp $p "${S}_rs_SBRef.nii.gz"
					fi
				done
				if [ ! -r "${S}_rs_SBRef.nii.gz" ]; then 
					echo "*** ERROR: could not locate ${S}_rs_SBRef.nii, or SBRef 6vols Pfile.nii"
				fi
			fi
			# -- RS-fMRI pePA Rev Run -->
			if [ ! -r "${S}_rs_pePA.nii.gz" ]; then
				for p in $(imglob $pfileDIR/*.nii.gz); do
					if [ -r "${p}.nii" ]; then
						pigz -v -b 512 -n -f -6 ${p}.nii
					fi
					H=$(3dinfo -history ${p}.nii.gz | grep '3dLRflip')
					if [ x"$H" == 'x' -a "$SubjAcqDate" -ge 20180101 ]; then
						echo "*** ERROR: RS-NII file=[$p] does NOT have required history containing '3dLRflip'"
						exit 2					
					elif [ x"$H" != 'x' -a "$SubjAcqDate" -lt 20180101 ]; then
						echo "*** ERROR: RS-NII file=[$p] incorrectly has a history containing '3dLRflip'"
						exit 2
					fi
					hdr=$(fslhd $p)
					NV=$(echo "$hdr" | grep '^dim4' | awk {'print $2'} | bc)
					TR=$(echo "$hdr" | grep '^pixdim4' | awk {'print $2'} | bc)
					#echo "-- found P=$p, tr=[$TR],nv=[$NV]"
					if [ "$NV" -eq 4 -a "$TR" == ".950000" ]; then
						echo "-- copying rev peAP RS-fMRI run = $p"
						$FSLDIR/bin/imcp $p "${S}_rs_pePA.nii.gz"
					fi
				done
				if [ ! -r "${S}_rs_pePA.nii.gz" ]; then 
					echo "*** ERROR: could not locate ${S}_rs_pePA.nii, or RevBlips 4vols Pfile.nii"
				fi
			fi
			# -- RS-fMRI peAP Main Run -->
			if [ ! -r "${S}_rs_peAP.nii.gz" ]; then
				for p in $(imglob $pfileDIR/*.nii.gz); do
					if [ -r "${p}.nii" ]; then
						pigz -v -b 512 -n -f -6 ${p}.nii
					fi
					H=$(3dinfo -history ${p}.nii.gz | grep '3dLRflip')
					if [ x"$H" == 'x' -a "$SubjAcqDate" -ge 20180101 ]; then
						echo "*** ERROR: RS-NII file=[$p] does NOT have required history containing '3dLRflip'"
						exit 2					
					elif [ x"$H" != 'x' -a "$SubjAcqDate" -lt 20180101 ]; then
						echo "*** ERROR: RS-NII file=[$p] incorrectly has a history containing '3dLRflip'"
						exit 2
					fi
					hdr=$(fslhd $p)
					NV=$(echo "$hdr" | grep '^dim4' | awk {'print $2'} | bc)
					TR=$(echo "$hdr" | grep '^pixdim4' | awk {'print $2'} | bc)
					#echo "-- found P=$p, tr=[$TR],nv=[$NV]"
					if [ "$NV" -eq 500 -a "$TR" == ".950000" ]; then
						echo "-- copying main peAP RS-fMRI run = $p"
						$FSLDIR/bin/imcp $p "${S}_rs_peAP.nii.gz"
					fi
				done
				if [ ! -r "${S}_rs_peAP.nii.gz" ]; then 
					echo "*** ERROR: could not locate ${S}_rs_peAP.nii, or Main 500 vols Pfile.nii"
				fi
			fi
		fi
	fi  ## if $SETUP_RS_PFILES ; then ...
	
	## re-orient NIFTI files to RPI
	if [ "$REQUIRE_REORIENT_CHECK" == "yes" ]; then
		for img in $(find . -maxdepth 1 -type f -name "*.nii.gz" -print); do
			o=$(3dinfo -orient ${img} | awk '{print $1}')
			if [ "$o" != "RPI" ] ; then
				nimg="$(basename $img .nii.gz)_${o}.nii.gz"
				mv -v $img $nimg
				3dresample -orient RPI -input $nimg -prefix $img
				rm -f $nimg
				echo " ++ reorienting $img to RPI"
			fi
		done
	fi
	
	echo "--------- finished `basename $0` on `pwd`, at: `date`"	
	ls -l 
done

