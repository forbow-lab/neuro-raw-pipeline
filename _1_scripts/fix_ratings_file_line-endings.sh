#!/bin/bash

cd /share/uher/FORBOW/rawdata/

for S in `/bin/ls -d ???_?_20??????`; do 
	echo "-----$S "
	ssid=`echo $S | awk -F_ '{print $1"_"$2}'`
	f="$S/${ssid}_qc_ratings"
	mv ${f}.csv ${f}_original.csv 
	tr '\r' '\n' <${f}_original.csv >${f}.csv
	numLines=$(cat ${f}.csv | wc -l)
	if [ "$numLines" -lt 3 ]; then
		rm ${f}.csv 
		mv ${f}_original.csv ${f}.csv 
	elif [ "$numLines" -eq 3 ]; then
		echo "" >>${f}.csv	##add extra line to file
		rm ${f}_original.csv
	else
		rm ${f}_original.csv
	fi
	cat ${f}.csv
done

exit 0
