#!/bin/bash

cd /shared/uher/FORBOW/rawdata/

for S in $(/bin/ls -d ???_?_201[8-9]????); do 
	for f in $(find $S/RS -type f -name "201?????_E*_P?????.7.nii.gz" -print); do 
		if [ "x$(3dinfo $f | grep '3dLRflip -X')" == "x" ] ; then 
			echo " $f  is not LRflipped"; 
		fi; 
	done; 
done

exit 0

