#!/bin/bash

export OMP_NUM_THREADS=4
export FSLDIR=/usr/local/fsl-5.0.10
source ${FSLDIR}/etc/fslconf/fsl.sh
export PATH="$FSLDIR/bin:$PATH" 
export FSLOUTPUTTYPE=NIFTI_GZ

REPORT_POSITIVE=false

Usage() {
	echo 
	echo "Usage: `basename $0`"
	echo 
	echo "Example: "
	echo "    cd /shared/uher/FORBOW/rawdata/"
	echo "    ./_1_scripts/compare_SSID_to_DICOM_id.sh"
	echo
	exit 1
}


#[ "$#" -lt 1 ] && Usage

SCRIPT=$(python -c "import os; print os.path.abspath('$0')")
SCRIPTSDIR=$(dirname $SCRIPT)
EDIR=/shared/uher/FORBOW/rawdata 

cd $EDIR/
SUBJECTS=$(/bin/ls -d ???_?_20??????)	#"$@"


for S in $SUBJECTS ; do 
	
	D=`find $S -maxdepth 1 -name "DICOMS" | head -1` 
	if [ ! -d "$D" ]; then 
		echo "$S  - no DICOMS directory found" 
		continue
	fi
	t1DIR=`find $D -maxdepth 1 -name "*T1w_BRAVO_sag*" | head -1`
	if [ ! -d "$t1DIR" ]; then 
		echo "$S -- no T1 directory found"
		continue
	fi
	im=`find $t1DIR -maxdepth 1 -name "IM-*-0001.dcm" | head -1`
	if [ ! -r "$im" ] ; then 
		echo "$S -- no T1 IM-*-0001.dcm found..."
		continue
	fi
	dicomMatch=false
	SSID=`echo $S | awk -F_ '{print $1"_"$2}'`;
	DID=`dicom_hdr $im | grep 'Patient ID' | awk -F// '{print $3}' | sed 's/ *//g'`
	if [ "$SSID" != "$DID" ]; then 
		echo "$SSID does NOT match DICOM-ID=$DID"
	else
		dicomMatch=true
	fi
	
	
	RS=`find $S -maxdepth 1 -name "RS" | head -1` 
	if [ ! -d "$RS" ]; then 
		echo "$S  - no RS directory" 
		continue
	fi
	isGZ=false
	P=`find $RS -maxdepth 1 -type f -name "P?????.7" -print0 | xargs -0 ls -Sr | head -1`
	if [ ! -r "$P" ]; then
		P=`find $RS -maxdepth 1 -type f -name "P?????.7.gz" -print0 | xargs -0 ls -Sr| head -1`
		if [ ! -r "$P" ]; then 
			echo "$S  - no Pfile found" 
			continue
		fi
		isGZ=true
		gunzip $P
		P="`dirname $P`/`basename $P .gz`"
	fi
	pfileMatch=false
	PID=`rdgehdr $P | grep 'Patient ID' | awk '{print $4}'`
	if [ "$SSID" != "$PID" ]; then 
		echo "$SSID does NOT match PFILE-ID=$PID"
	else
		pfileMatch=true
	fi
	if $isGZ ; then
		pigz $P
	fi
	#find $RS -type f -name "P?????.7" -print0 | xargs pigz
	if $REPORT_POSITIVE && $dicomMatch && $pfileMatch ; then
		echo "$SSID matches both DICOM=$DID and PFILE=$PID"
	fi
done


exit 0
