#!/bin/bash

# To run: 
#  $ cd /shared/uher/FORBOW/rawdata/
#  $ D=$(date +%Y%m%d)
#  $ ./_1_scripts/FORBOW_confirm_complete_bids.sh 2>./_3_docs/${D}_confirm_bids_errors.txt 1>./_3_docs/${D}_confirm_bids_ready.txt
#

if [ "`uname -n`" == "Aoraki.local" ]; then
	PROJECT_DIR="$HOME/work/uher/FORBOW"
else
	PROJECT_DIR="/shared/uher/FORBOW"
fi

cd $PROJECT_DIR/rawdata/
for D in `/bin/ls -d ???_?_20??????`; do
	S=${D:0:3}
	E=${D:4:1}
	BD="$D/BIDS_raw/sub-${S}/ses-${E}"
	AD="$BD/anat"
	DD="$BD/dwi"
	FMD="$BD/fmap"
	FUD="$BD/func"
	pf="sub-${S}_ses-${E}"
	goodBIDS="yes"
	for f in $AD/${pf}_T1w.nii.gz $AD/${pf}_T1w.json $AD/${pf}_mod-T1w_defacemask.nii.gz $AD/${pf}_T2w.nii.gz $AD/${pf}_T2w.json $AD/${pf}_mod-T2w_defacemask.nii.gz ; do
		if [ ! -f "$f" ]; then
			echo " - missing $f" >&2
			goodBIDS="no"
		fi
	done
	for pe in AP PA ; do
		for f in $DD/${pf}_dir-${pe}_dwi.nii.gz $DD/${pf}_dir-${pe}_dwi.json $DD/${pf}_dir-${pe}_dwi.bval $DD/${pf}_dir-${pe}_dwi.bvec ; do
			if [ ! -f "$f" ]; then
				echo " - missing $f" >&2
				goodBIDS="no"
			fi
		done
	done
	for f in $FMD/${pf}_dir-AP_epi.nii.gz $FMD/${pf}_dir-AP_epi.json $FMD/${pf}_dir-PA_epi.nii.gz $FMD/${pf}_dir-PA_epi.json ; do
		if [ ! -f "$f" ]; then
			echo " - missing $f" >&2
			goodBIDS="no"
		fi
	done
	for f in $FUD/${pf}_task-rest_bold.nii.gz $FUD/${pf}_task-rest_bold.json $FUD/${pf}_task-rest_sbref.nii.gz $FUD/${pf}_task-rest_sbref.json ; do
		if [ ! -f "$f" ]; then
			echo " - missing $f" >&2
			goodBIDS="no"
		fi
	done
	if [ "$goodBIDS" == "yes" ]; then
		echo " + bids complete = $BD" >&1
	fi
done
exit 0



BIDS_raw/README
BIDS_raw/dataset_description.json
BIDS_raw/task-rest_bold.json

BIDS_raw/sub-131/sub-131_sessions.tsv

BIDS_raw/sub-131/ses-D/anat/
sub-131_ses-D_T1w.json
sub-131_ses-D_T1w.nii.gz
sub-131_ses-D_T2w.json
sub-131_ses-D_T2w.nii.gz
sub-131_ses-D_mod-T1w_defacemask.nii.gz
sub-131_ses-D_mod-T2w_defacemask.nii.gz

BIDS_raw/sub-131/ses-D/dwi/
sub-131_ses-D_dir-AP_dwi.bval
sub-131_ses-D_dir-AP_dwi.bvec
sub-131_ses-D_dir-AP_dwi.json
sub-131_ses-D_dir-AP_dwi.nii.gz
sub-131_ses-D_dir-PA_dwi.bval
sub-131_ses-D_dir-PA_dwi.bvec
sub-131_ses-D_dir-PA_dwi.json
sub-131_ses-D_dir-PA_dwi.nii.gz

BIDS_raw/sub-131/ses-D/fmap/
sub-131_ses-D_dir-AP_epi.json
sub-131_ses-D_dir-AP_epi.nii.gz
sub-131_ses-D_dir-PA_epi.json
sub-131_ses-D_dir-PA_epi.nii.gz

BIDS_raw/sub-131/ses-D/func/
sub-131_ses-D_task-rest_bold.json
sub-131_ses-D_task-rest_bold.nii.gz
sub-131_ses-D_task-rest_sbref.json
sub-131_ses-D_task-rest_sbref.nii.gz



# for D in `ls -d ???_?_20??????` ; do S=${D:0:5}; if [ -r "$D/NIFTIS/${S}_t1w_bravo_orig_defaced.nii.gz" -a -r "$D/NIFTIS/${S}_t2prep_flair_promo_orig_defaced.nii.gz" ]; then echo " + $D has both T1w_defaced and T2Prep_defaced"; fi;  done

if [ "`uname -n`" == "Aoraki.local" ]; then
	PROJECT_DIR="$HOME/work/uher/FORBOW"
else
	PROJECT_DIR="/shared/uher/FORBOW"
fi

cd $PROJECT_DIR/rawdata/
for D in $(/bin/ls -d ???_?_20??????) ; do 
	S=${D:0:5};
	hasT1w=yes; hasT2Flair=yes; hasDWI=yes; hasRS=yes
	if [ ! -r "$D/NIFTIS/${S}_t1w_bravo_orig.nii.gz" ]; then 
		echo " --- $S missing t1w_bravo_orig.nii" >&2
		hasT1w=no
	fi
	if [ ! -r "$D/NIFTIS/${S}_t2prep_flair_promo_orig.nii.gz" ]; then
		echo " --- $S missing t2prep_flair_promo_orig.nii" >&2
		hasT2Flair=no
	fi
	if [ ! -r "$D/NIFTIS/${S}_dwi_peAP.nii.gz" -o ! -r "$D/NIFTIS/${S}_dwi_pePA.nii.gz" ]; then
		echo " --- $S missing dwi_peAP.nii OR dwi_pePA.nii" >&2
		hasDWI=no
	fi
	if [ ! -r "$D/NIFTIS/${S}_rs_peAP.nii.gz" -o ! -r "$D/NIFTIS/${S}_rs_pePA.nii.gz" ]; then
		echo " --- $S missing rs_peAP.nii OR rs_pePA.nii" >&2
		hasRS=no
	fi
	if [ "$hasT1w" == "yes" -a "$hasT2Flair" == "yes" -a "$hasDWI" == "yes" -a "$hasRS" == "yes" ] ; then
		echo " + $S has all required sequences"
	fi
done

exit 0