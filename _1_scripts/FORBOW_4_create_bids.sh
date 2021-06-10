#!/bin/bash
#-----------------------------------------------------------------------------------------
# 2020/12/09 - Carl Helmick
#	* need to create SESS variable separate from SSID? 
#	* copy T1.json file header to get scanner-details, etc.
#	* add all extra RS sequence information necessary for BIDS-compliance 
#	* match PhaseEncodingDirection="j" of DWI-sequence from DICOMS dcm2nix
#     and using visual comparison between DWI and RS nifti files.
#
#  ASSUMPTIONS:
#	1) GE collects axial scans as R,P,I as 0,0,0. Therefore:
#		PEpolarity of 1  = P->A = j
#		PEpolarity of -1 = A->P = j-
#     Q) Does this actual direction label matter for correction of BIDS, or only consistent usage within a study?
#
# MODIFICATIONS:
#   2021/02/24, Carl:
#      + updating project paths
#----------------------------------------------------------------------------------------

Usage(){
	echo "Usage: `basename $0` <ssid>"
	echo
	echo "Example: `basename $0` 005_C"
	echo
	exit 1
}

SCRIPT=$(python -c "from os.path import abspath; print(abspath('$0'))")
SCRIPTSDIR=$(dirname $SCRIPT)

## MAIN PATH to PROJECT ###########
if [ "$HOSTNAME" == "Aoraki.local" ]; then
	PROJECT_DIR="$HOME/work/uher/FORBOW"
else
	PROJECT_DIR="/shared/uher/FORBOW"
fi
RAW_DATA_DIR="${PROJECT_DIR}/rawdata"
SOFTWARE_DIR="${RAW_DATA_DIR}/_0_software"
BIDS_TEMPLATE_DIR="${SCRIPTSDIR}/BIDS_Templates"
DWI_pePA_NUMVOLS="33"
DWI_peAP_NUMVOLS="8"
FORCE_OVERWRITE="yes"


[ "$#" -lt 1 ] && Usage


SUBJECTS=""
for Subj in $@ ; do
	
	S=$(basename $Subj)
	if [ "${#S}" -eq 14 ]; then
		SSID="${S:0:3}"
		SESS="${S:4:1}"
		SubjRawDIR="${RAW_DATA_DIR}/$S"
		S="${SSID}_${SESS}"
		if [ ! -d "$SubjRawDIR" ]; then
			echo "*** ERROR: cannot find rawdata directory for subject=[$S] in ${RAW_DATA_DIR}/"
			continue
		fi
	elif [ "${#S}" -eq 5 ]; then
		SubjRawDIR=$(find ${RAW_DATA_DIR} -maxdepth 1 -type d -name "${S}_20??????" -print)
		if [ ! -d "$SubjRawDIR" ]; then
			echo "*** ERROR: cannot find rawdata directory for subject=[$S] in ${RAW_DATA_DIR}/"
			continue
		fi
		SSID="${S:0:3}"
		SESS="${S:4:1}"
	else
		echo " * ERROR: cannot determine subject from input=[$Subj]"
		continue
	fi
	SDIR_BN=$(basename $SubjRawDIR)
	if [ "${#SDIR_BN}" -eq 14 ]; then
		SubjAcqDate=${SDIR_BN:6:8}
	else
		echo "*** ERROR: failed to determine AcqDate from SubjRawDIR=${SubjRawDIR}..."
		continue
	fi
	if [ ! -d "$SubjRawDIR" ]; then
		echo "*** ERROR: cannot find rawdata directory for subject=[$S] in ${RAW_DATA_DIR}/"
		continue
	fi
	SubjNiiDIR="$SubjRawDIR/NIFTIS"
	if [ ! -d "$SubjNiiDIR" ]; then
		if [ ! -d "$SubjNiiDIR" ]; then
			echo "*** ERROR: could not find directory = $SubjNiiDIR/"
			continue
		fi
	fi
	SubjDicomDIR="$SubjRawDIR/DICOMS"
	if [ ! -d "$SubjDicomDIR" ]; then
		echo "*** ERROR: could not find directory = $SubjDicomDIR/"
		continue
	fi
	SubjPfileDIR="$SubjRawDIR/RS"
	if [ ! -d "$SubjPfileDIR" ]; then
		echo "*** ERROR: could not find directory = $SubjPfileDIR/"
		continue
	fi


	echo " ---> running $(basename $0) on ${SubjRawDIR}, starting `date`"	
	
	SubjBIDSdir="${SubjRawDIR}/BIDS_raw"
	if [ "$FORCE_OVERWRITE" == "yes" -a -d "$SubjBIDSdir" ]; then
		echo " * FORCE_OVERWRITE enabled, removing bidsdir = $SubjBIDSdir"
		rm -rf ${SubjBIDSdir}/
	fi
	
	
	SubjBIDS="${SubjRawDIR}/BIDS_raw/sub-${SSID}"
	SubjSessBIDS="${SubjBIDS}/ses-${SESS}"
	mkdir -p ${SubjSessBIDS}/anat/
	newSubj="no"
	## convert T1w ----------------------------------------------------------------------
	if [ ! -r "${SubjSessBIDS}/anat/sub-${SSID}_ses-${SESS}_T1w.nii.gz" ]; then 
	
		## FORBOW has ${S}_t1w_bravo_orig.nii.gz
		if [ ! -r "$SubjNiiDIR/${S}_t1w_bravo_orig.nii.gz" ]; then
			echo "*** ERROR: could not locate subject file = $SubjNiiDIR/${S}_t1w_bravo_orig.nii.gz"
			continue
		fi
		cp -pv $SubjNiiDIR/${S}_t1w_bravo_orig.nii.gz ${SubjSessBIDS}/anat/sub-${SSID}_ses-${SESS}_T1w.nii.gz
		cp -pv $SubjNiiDIR/${S}_t1w_bravo_orig.json ${SubjSessBIDS}/anat/sub-${SSID}_ses-${SESS}_T1w.json
	fi
	## https://bids-specification.readthedocs.io/en/stable/99-appendices/09-entities.html#mod
	if [ ! -r "${SubjSessBIDS}/anat/sub-${SSID}_ses-${SESS}_mod-T1w_defacemask.nii.gz" ]; then 
		if [ ! -r "$SubjNiiDIR/${S}_t1w_bravo_orig_facemask.nii.gz" ]; then
			echo "*** ERROR: could not locate subject file = $SubjNiiDIR/${S}_t1w_bravo_orig_facemask.nii.gz"
			continue
		fi
		cp -pv $SubjNiiDIR/${S}_t1w_bravo_orig_facemask.nii.gz ${SubjSessBIDS}/anat/sub-${SSID}_ses-${SESS}_mod-T1w_defacemask.nii.gz
	fi
	
	## convert T2w ---------------------------------------------------------------------
	if [ ! -r "${SubjSessBIDS}/anat/sub-${SSID}_ses-${SESS}_T2w.nii.gz" ]; then
		## FORBOW has ${S}_t2prep_flair_promo_orig.nii.gz
		if [ ! -r "$SubjNiiDIR/${S}_t2prep_flair_promo_orig.nii.gz" ]; then
			echo "*** ERROR: could not locate subject file = $SubjNiiDIR/${S}_t2prep_flair_promo_orig.nii.gz"
			continue
		fi
		cp -pv $SubjNiiDIR/${S}_t2prep_flair_promo_orig.nii.gz ${SubjSessBIDS}/anat/sub-${SSID}_ses-${SESS}_T2w.nii.gz
		cp -pv $SubjNiiDIR/${S}_t2prep_flair_promo_orig.json ${SubjSessBIDS}/anat/sub-${SSID}_ses-${SESS}_T2w.json
	fi
	## https://bids-specification.readthedocs.io/en/stable/99-appendices/09-entities.html#mod
	if [ ! -r "${SubjSessBIDS}/anat/sub-${SSID}_ses-${SESS}_mod-T2w_defacemask.nii.gz" ]; then 
		if [ ! -r "$SubjNiiDIR/${S}_t2prep_flair_promo_orig_facemask.nii.gz" ]; then
			echo "*** ERROR: could not locate subject file = $SubjNiiDIR/${S}_t2prep_flair_promo_orig_facemask.nii.gz"
			continue
		fi
		cp -pv $SubjNiiDIR/${S}_t2prep_flair_promo_orig_facemask.nii.gz ${SubjSessBIDS}/anat/sub-${SSID}_ses-${SESS}_mod-T2w_defacemask.nii.gz
	fi
	
	
	## create subject-session file if none exists
	FirstDicom=$(find $SubjDicomDIR -maxdepth 3 -type f -iname "IM-000*-0001.dcm" -print | head -1)
	if [ ! -r "$FirstDicom" ]; then
		echo "*** ERROR: could not locate an example DICOM file to determine scan date/time..."
	else
		sessFile=${SubjBIDS}/sub-${SSID}_sessions.tsv
		SubjAcqTime=$(dicom_hdr $FirstDicom | grep 'ID Study Time' | awk -F// '{print $3}')
		SubjAcqDate=$(dicom_hdr $FirstDicom | grep 'ID Study Date' | awk -F// '{print $3}')
		acq_time="${SubjAcqDate:0:4}-${SubjAcqDate:4:2}-${SubjAcqDate:6:2}T${SubjAcqTime:0:2}:${SubjAcqTime:2:2}:${SubjAcqTime:4:2}"
		if [ ! -r "${sessFile}" ]; then
			echo -e "session_id\tacq_time" >${sessFile}
			echo -e "ses-${SESS}\t${acq_time}" >>${sessFile}
		else
			sessListed=$(cat -v ${sessFile} | grep "ses-$SESS")
			if [ "$sessListed" == "" ]; then
				echo -e "ses-${SESS}\t${acq_time}" >>${sessFile}
			fi
		fi
	fi
	
	## add task-bold.json, README, and
	cp -vpf $BIDS_TEMPLATE_DIR/dataset_description.json  $SubjRawDIR/BIDS_raw/
	cp -vpf $BIDS_TEMPLATE_DIR/task-rest_bold.json  $SubjRawDIR/BIDS_raw/
	cp -vpf $BIDS_TEMPLATE_DIR/README  $SubjRawDIR/BIDS_raw/
	
	
	#####  DWI  ------------------------------------------------------------------------
	mkdir -p ${SubjSessBIDS}/dwi/
	## convert DWI_MAIN
	dwi="sub-${SSID}_ses-${SESS}_dir-PA_dwi"
	if [ ! -r "${SubjSessBIDS}/dwi/${dwi}.nii.gz" ]; then
		if [ ! -r "$SubjNiiDIR/${S}_dwi_pePA.nii.gz" ]; then
			echo "*** ERROR: could not locate subject file = $SubjNiiDIR/${S}_dwi_pePA.nii.gz"
			continue
		fi
		cp -pv $SubjNiiDIR/${S}_dwi_pePA.nii.gz ${SubjSessBIDS}/dwi/${dwi}.nii.gz
		cp -pv $SubjNiiDIR/${S}_dwi_pePA.json ${SubjSessBIDS}/dwi/${dwi}.json
		cp -pv $SubjNiiDIR/${S}_dwi_pePA.bvec ${SubjSessBIDS}/dwi/${dwi}.bvec
		cp -pv $SubjNiiDIR/${S}_dwi_pePA.bval ${SubjSessBIDS}/dwi/${dwi}.bval
	fi
	## convert DWI_REV
	dwi="sub-${SSID}_ses-${SESS}_dir-AP_dwi"
	if [ ! -r "${SubjSessBIDS}/dwi/${dwi}.nii.gz" ]; then
		if [ ! -r "$SubjNiiDIR/${S}_dwi_peAP.nii.gz" ]; then
			echo "*** ERROR: could not locate subject file = $SubjNiiDIR/${S}_dwi_peAP.nii.gz"
			continue
		fi
		cp -pv $SubjNiiDIR/${S}_dwi_peAP.nii.gz ${SubjSessBIDS}/dwi/${dwi}.nii.gz
		cp -pv $SubjNiiDIR/${S}_dwi_peAP.json ${SubjSessBIDS}/dwi/${dwi}.json
		cp -pv $SubjNiiDIR/${S}_dwi_peAP.bvec ${SubjSessBIDS}/dwi/${dwi}.bvec
		cp -pv $SubjNiiDIR/${S}_dwi_peAP.bval ${SubjSessBIDS}/dwi/${dwi}.bval
	fi
	
	
	
	#####  func && fmap -------------------------------------------------------------------------
	mkdir -p ${SubjSessBIDS}/func/ ${SubjSessBIDS}/fmap/ && cd ${SubjSessBIDS}/func/
	## from NSEPP-CIHR-CB053
	#	
	#  dwi_pePA.json: "PhaseEncodingDirection":"j" ,"PhaseEncodingPolarityGE":"Flipped"==DTI_30_DIR+3B0==P->A
	#  dwi_peAP.json: "PhaseEncodingDirection":"j-","PhaseEncodingPolarityGE":"Unflipped"==DTI_B0map__Rev==A->P
	#
	#  $ for f in *_hdr.txt ; do echo "$f: `cat -v $f | grep 'rhdacqctrl'`"; done
	#  >20190904-0612_E14874_P16384.7_hdr.txt: ...rhdacqctrl bit mask: 8	rs_SBref
	#  >20190904-0612_E14874_P16896.7_hdr.txt: ...rhdacqctrl bit mask: 12	rs_peAP
	#  >20190904-0612_E14874_P17408.7_hdr.txt: ...rhdacqctrl bit mask: 8	rs_pePA
	#
	# TotalNumberVols: 
	# typical mux_epi_recon ignores first 6vols when doing recon to allow for field to reach steady-state.
	# This does mean that Pfiles will report 6extra volumes (nt) than output NIFTI file.
	#
	## RS_pePA - 500vols_fwd
	## this is the main rs run with vols=500, TR=0.950,PE=PA, 3mm-isotropic, MB=3, matrix=72x72x71 
	nVolsReq="500"
	TRreq="0.950000"
	PEdir="PA"
	func="${SubjSessBIDS}/func/sub-${SSID}_ses-${SESS}_task-rest_bold"
	if [ ! -r "${func}.json" ]; then		
		echo " + searching for NIFTI file with ${nVolsReq}vols with TR=${TRreq}sec "
		for pfile in $(find $SubjPfileDIR -maxdepth 1 -type f -name "20*_E*_P*.nii.gz" -print); do
			TR=$(fslval $pfile pixdim4 | awk '{x=$1; printf "%.06f",x}')
			nVols=$(fslval $pfile dim4 | awk '{x=$1; printf "%d",x}')
			if [ "$TR" == "$TRreq" -a "$nVols" == "$nVolsReq" ]; then
				echo " ++ found ${nVols}vol run with TR=${TR}: $pfile"
				cp -pv $pfile ${func}.nii.gz
				tmpJson="${SubjSessBIDS}/dwi/sub-${SSID}_ses-${SESS}_dir-${PEdir}_dwi.json"
				echo " + using as template json=${tmpJson}"
				if [ ! -r "$tmpJson" ]; then
					echo "*** ERROR: cannot locate json-file = $tmpJson"
					exit 7
				fi
				pfileHDR="${pfile:0:$((${#pfile}-7))}_hdr.txt"
				if [ ! -r "$pfileHDR" ]; then
					echo "*** ERROR: cannot locate pfile_hdr = $pfileHDR"
					exit 7
				fi
				pfileMAT="${pfile:0:$((${#pfile}-7))}.mat"
				if [ ! -r "$pfileMAT" ]; then
					echo "*** ERROR: cannot locate pfile_mat = $pfileMAT"
					exit 7
				fi
				outJsonFile="${func}.json"
				echo " + running cmd: ${SCRIPTSDIR}/create_json_for_Pfile.py -j $tmpJson -p $pfileHDR -m $pfileMAT -o $outJsonFile"
				${SCRIPTSDIR}/bids_create_json_for_Pfile.py -j $tmpJson -p $pfileHDR -m $pfileMAT -o $outJsonFile
				if [ "$?" -ne 0 ]; then
					echo "*** ERROR: was unable to create json-file = $outJsonFile"
					exit 7
				else
					## now also create the fmap
					mkdir -p ${SubjSessBIDS}/fmap/
					fmap="${SubjSessBIDS}/fmap/sub-${SSID}_ses-${SESS}_dir-${PEdir}_epi"
					echo " ++ mcflirt -meanvol for motion-correction, with fslmaths -Tmean, creating: ${fmap}.nii"
					mcflirt -in ${func}.nii -meanvol -spline_final -o ${fmap}
					fslmaths ${fmap}.nii -thr 0 -Tmean ${fmap}
					rm -f ${fmap}_mean_reg.nii.gz
					cat -v $outJsonFile | grep -e "{" -e "PhaseEncodingPolarityGE" -e "EffectiveEchoSpacing" -e "TotalReadoutTime" -e "PhaseEncodingDirection" -e "MultibandAccelerationFactor" >${fmap}.json
					echo -e "    \"IntendedFor\": \"ses-${SESS}/func/sub-${SSID}_ses-${SESS}_task-rest_bold.nii.gz\"\n}\n" >>${fmap}.json
				fi				
			fi
		done
	fi
	
	## RS_peAP - 4vols_rev
	## this is the rev rs run with vols=4, TR=0.950, PE=AP, 3mm-isotropic, MB=3, matrix=72x72x71 
	nVolsReq="4"
	TRreq="0.950000"
	PEdir="AP"
	ACQmux="mb3"
	fmap="${SubjSessBIDS}/fmap/sub-${SSID}_ses-${SESS}_dir-${PEdir}_epi"
	if [ ! -r "${fmap}.json" ]; then		
		echo " + searching for NIFTI file with ${nVolsReq}vols with TR=${TRreq}sec "
		for pfile in $(find $SubjPfileDIR -maxdepth 1 -type f -name "20*_E*_P*.nii.gz" -print); do
			TR=$(fslval $pfile pixdim4 | awk '{x=$1; printf "%.06f",x}')
			nVols=$(fslval $pfile dim4 | awk '{x=$1; printf "%d",x}')
			if [ "$TR" == "$TRreq" -a "$nVols" == "$nVolsReq" ]; then
				echo " ++ found ${nVols}vol run with TR=${TR}: $pfile"
				cp -pv $pfile ${fmap}_raw.nii.gz
				tmpJson=${SubjSessBIDS}/dwi/sub-${SSID}_ses-${SESS}_dir-${PEdir}_dwi.json
				pfileHDR="${pfile:0:$((${#pfile}-7))}_hdr.txt"
				echo " + using as template json=${tmpJson}"
				if [ ! -r "$pfileHDR" ]; then
					echo "*** ERROR: cannot locate pfile_hdr = $pfileHDR"
					exit 7
				fi
				pfileMAT="${pfile:0:$((${#pfile}-7))}.mat"
				if [ ! -r "$pfileMAT" ]; then
					echo "*** ERROR: cannot locate pfile_mat = $pfileMAT"
					exit 7
				fi
				mkdir -p ${SubjSessBIDS}/fmap/
				${SCRIPTSDIR}/bids_create_json_for_Pfile.py -j $tmpJson -p $pfileHDR -m $pfileMAT -o ${fmap}_tmp.json
				if [ "$?" -ne 0 ]; then
					echo "*** ERROR: was unable to create json-file = ${fmap}_tmp.json"
					exit 7
				else
					## now also create the fmap
					cat -v ${fmap}_tmp.json | grep -e "{" -e "PhaseEncodingPolarityGE" -e "EffectiveEchoSpacing" -e "TotalReadoutTime" -e "PhaseEncodingDirection" -e "MultibandAccelerationFactor" >${fmap}.json
					echo -e "    \"IntendedFor\": \"ses-${SESS}/func/sub-${SSID}_ses-${SESS}_task-rest_bold.nii.gz\"\n}\n" >>${fmap}.json
					rm -f ${fmap}_tmp.json
					echo " ++ mcflirt -meanvol for motion-correction, with fslmaths -Tmean, creating: ${fmap}.nii"
					mcflirt -in ${fmap}_raw.nii.gz -meanvol -spline_final -o ${fmap}
					fslmaths ${fmap}.nii -thr 0 -Tmean ${fmap}
					rm -f ${fmap}_mean_reg.nii.gz ${fmap}_raw.nii.gz
				fi
			fi
		done
	fi
	
	## RS_SBref - 6vols_fwd
	## this is the SBref rs run with vols=6, TR=[2.8-3.0], PE=PA, 3mm-isotropic, MB=1, matrix=72x72x71 
	nVolsReq="6"
	TRreq="2.800000"
	PEdir="PA"
	ACQmux="sb1"
	func="${SubjSessBIDS}/func/sub-${SSID}_ses-${SESS}_task-rest_sbref"
	if [ ! -r "${func}.json" ]; then		
		echo " + searching for NIFTI file with ${nVolsReq}vols with TR=${TRreq}sec "
		for pfile in $(find $SubjPfileDIR -maxdepth 1 -type f -name "20*_E*_P*.nii.gz" -print); do
			TR=$(fslval $pfile pixdim4 | awk '{x=$1; printf "%.06f",x}')
			nVols=$(fslval $pfile dim4 | awk '{x=$1; printf "%d",x}')
			if [ "$TR" == "$TRreq" -a "$nVols" == "$nVolsReq" ]; then
				echo " ++ found ${nVols}vol run with TR=${TR}: $pfile"
				cp -pv $pfile ${func}_raw.nii.gz
				echo " ++ mcflirt -meanvol for motion-correction, with fslmaths -Tmean, creating: ${fmap}.nii"
				mcflirt -in ${func}_raw.nii -meanvol -spline_final -o ${func}
				fslmaths ${func}.nii -thr 0 -Tmean ${func}
				rm -f ${func}_mean_reg.nii.gz ${func}_raw.nii.gz
				tmpJson=${SubjSessBIDS}/dwi/sub-${SSID}_ses-${SESS}_dir-${PEdir}_dwi.json
				pfileHDR="${pfile:0:$((${#pfile}-7))}_hdr.txt"
				echo " + using as template json=${tmpJson}"
				if [ ! -r "$pfileHDR" ]; then
					echo "*** ERROR: cannot locate pfile_hdr = $pfileHDR"
					exit 7
				fi
				pfileMAT="${pfile:0:$((${#pfile}-7))}.mat"
				if [ ! -r "$pfileMAT" ]; then
					echo "*** ERROR: cannot locate pfile_mat = $pfileMAT"
					exit 7
				fi
				outJsonFile=${func}.json
				${SCRIPTSDIR}/bids_create_json_for_Pfile.py -j $tmpJson -p $pfileHDR -m $pfileMAT -o $outJsonFile
				if [ "$?" -ne 0 ]; then
					echo "*** ERROR: was unable to create json-file = ${outJsonFile}"
					exit 7
				fi
			fi
		done
	fi
	
	
	echo " ---> completed running $(basename $0) on ${SubjRawDIR}, `date`"	
		
done


#cd ${BIDS_RAW_DIR}/
#echo " ++ running subject-level MRIQC on ${S}"
#docker run -it --rm -v ${BIDS_RAW_DIR}:/data:ro -v ${BIDS_RAW_DIR}/mriqc_results:/out poldracklab/mriqc:latest /data /out participant --no-sub --participant_label ${SUBJECTS}


exit 0
