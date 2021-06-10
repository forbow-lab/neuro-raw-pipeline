#!/bin/bash

Usage() {
	echo
	echo "Usage: `basename $0` <ssid_sess>"
	echo
	echo "Example: `basename $0` 018_C 123_B"
	echo
	exit 1
}

SCRIPT=$(python -c "from os.path import abspath; print(abspath('$0'))")
SCRIPTSDIR=$(dirname $SCRIPT)
if [[ -z "${RAW_DATA_DIR}" ]]; then
	source "$SCRIPTSDIR/FORBOW_SetupEnvironment.sh"
fi

run_full_pipeline() {
	echo "--------------------------------------------------------------------"
	echo "`date`: starting `basename $0` on subject=$1"
	$SCRIPTSDIR/FORBOW_1_run_mux2nii.sh $1
	$SCRIPTSDIR/FORBOW_2_create_niftis.sh $1
	$SCRIPTSDIR/FORBOW_3_deface_and_QC.sh $1
	$SCRIPTSDIR/FORBOW_4_create_bids.sh $1
	echo "`date`: finished `basename $0` on subject=$1" 
	echo "--------------------------------------------------------------------"
}

[ "$#" -lt 1 ] && Usage
SUBJECTS="$@"
for Subj in ${SUBJECTS} ; do 
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
	if [ ! -d "$SDIR" ]; then
		echo "*** ERROR: cannot find rawdata directory for $S in ${RAW_DATA_DIR}/"
		continue
	fi
	echo " -- starting `basename $0` on $SDIR/, `date`"
	mkdir -p ${SDIR}/logs/
	LOG="${SDIR}/logs/${S}_rawdata_pipeline.txt"
	run_full_pipeline ${S} >>$LOG 2>&1
	echo " -- completed `basename $0` on $SDIR/, `date`"
done

exit 0
