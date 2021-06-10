#!/bin/bash


function report_missing(){

BadSubjs=""
for S in $@ ; do 
	t2FP=$(find $S/DICOMS -maxdepth 1 -type d -name "T2PREP_FLAIR_PROMO_*" -print | tail -1); 
	t1wO=$(find $S/DICOMS -maxdepth 1 -type d -name "ORIG_T1w_BRAVO_sag_*" -print | tail -1);  
	if [ -d "$t2FP" -a -d "$t1wO" ]; then 
		#echo "$S, $t2FP, $t1wO";
		continue 
	elif [ -d "$t2FP" -a ! -d "t1wO" ]; then
		echo "--- missing $S ORIG_T1w_BRAVO_sag"
		BadSubjs="${BadSubjs} $S"
	elif [ ! -d "$t2FP" -a -d "t1wO" ]; then
                echo "--- missing $S T2PREP_FLAIR_PROMO"
                BadSubjs="${BadSubjs} $S"
	else 
		echo "--- missing $S T2PREP_FLAIR_PROMO or ORIG_T1w_BRAVO_sag"; 
		BadSubjs="${BadSubjs} $S"
	fi
done

echo "BadSubjs = ${BadSubjs}"

}



cd /shared/uher/FORBOW/rawdata/

echo; echo " -- 2016 datasets ------------------------------------------->"
report_missing $(/bin/ls -d ???_?_2016????)

echo; echo " -- 2017 datasets ------------------------------------------->"
report_missing $(/bin/ls -d ???_?_2017????)

echo; echo " -- 2018 datasets ------------------------------------------->"
report_missing $(/bin/ls -d ???_?_2018????)

echo; echo " -- 2019 datasets ------------------------------------------->"
report_missing $(/bin/ls -d ???_?_2019????)



exit 0
