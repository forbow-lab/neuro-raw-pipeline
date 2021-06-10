#!/bin/bash


cd /shared/uher/FORBOW/rawdata/

for SDIR in `ls -d 00?_[A]_201[6-7]????`; do
	S=${SDIR:0:5}
	SCAN_YR=${SDIR:6:4}
	if [ ! -d "$SDIR/RS" ]; then
		echo "*** ERROR: cannot find rawdata/SSID/RS/ folder for $S"; 
		continue
	else
		echo "--------------- $SDIR/RS/"; 
		/bin/ls -l $SDIR/RS/; 
		for P in `/bin/ls -f $SDIR/RS/*.nii.gz`; do 
			LR=$(3dinfo -history $P | grep '3dLRflip -X');
			RPI=$(3dinfo -history $P | grep '3dresample -orient RPI')
			echo  "$P: rpi=${RPI},  lrflip=${LR}";
		done
	fi
	if [ ! -d "$SDIR/NIFTIS" ]; then 
		echo " - missing $SDIR/NIFTIS/" 
	fi
done
