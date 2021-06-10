#!/bin/bash

EXPDIR="/shared/uher/FORBOW"
BIDSDIR="$EXPDIR/mriqc/3T_BIDS"
RAWDIR="$EXPDIR/rawdata"

cd $RAWDIR/

    #for SD in `/bin/ls -d 042_[A,B,C,D,E,F]_*`; do 
for SD in `/bin/ls -d [0-1]??_[A,B,C,D,E,F]_*`; do 
    S=${SD:0:5}
    SSID=${S:0:3}
    SESS=${S:4:1}
    qcDIR=$BIDSDIR/sub-${SSID}/ses-${SESS}/anat

    if [ ! -f "$qcDIR/sub-${SSID}_ses-${SESS}_T1w.nii.gz" ]; then
        seqDIR=$(find $SD/DICOMS -maxdepth 1 -type d -name "T1w_BRAVO_sag_*" -print | tail -1)
        if [ -d "$seqDIR" ]; then   
            mkdir -p $qcDIR/
            echo "$SD --> $S --> $SSID-$SESS, T1wDIR=$seqDIR/"
            dcm2niix -b y -x n -z y -o $qcDIR/ -f "sub-${SSID}_ses-${SESS}_T1w" $seqDIR/
        fi
    fi

    if [ ! -f "$qcDIR/sub-${SSID}_ses-${SESS}_T2w.nii.gz" ]; then
        seqDIR=$(find $SD/DICOMS -maxdepth 1 -type d -name "T2w_CUBE_sag_*" -print | tail -1)
        if [ -d "$seqDIR" ]; then   
            mkdir -p $qcDIR/
            echo "$SD --> $S --> $SSID-$SESS, T2wDIR=$seqDIR/"
            dcm2niix -b y -x n -z y -o $qcDIR/ -f "sub-${SSID}_ses-${SESS}_T2w" $seqDIR/
        fi
    fi
     
    if [ ! -f "$qcDIR/sub-${SSID}_ses-${SESS}_FLAIR.nii.gz" ]; then
        seqDIR=$(find $SD/DICOMS -maxdepth 1 -type d -name "T2PREP_FLAIR_*" -print | tail -1)
        if [ -d "$seqDIR" ]; then   
            mkdir -p $qcDIR/
            echo "$SD --> $S --> $SSID-$SESS, FLAIRwDIR=$seqDIR/"
            dcm2niix -b y -x n -z y -o $qcDIR/ -f "sub-${SSID}_ses-${SESS}_FLAIR" $seqDIR/
        fi
    fi

done

exit 0

