#!/bin/bash
# 
# biotic_organize_pfiles.sh
#
# REQUIREMENTS: 
#	Linux OS with 'rdgehdr', 'pigz', 'pbzip2' and 'rename.ul'
#
# This script is meant to run on a directory containing pfiles, and will attempt three things:
# 1) Using rdgehdr, dumps Pfile header into text file = ${pfile}_hdr.txt
# 2) Renames each Pfile and all associated files (_vrgf.dat, _ref.h5, etc) to include a prefix of Exam date-time and exam number.
# 3) Compresses Pfile using pbzip2 (Pfiles will have .bz2 extension), achieving ~75% compression typically.
# -- Carl Helmick, 20180423
#
# Modifications:
#	20181018: 
# 		- autodetect and set NCORES for compressing/decompressing  = NUM_PHYSICAL_CORES
#		- replaced rename command with rename.ul
#-----------------------------------------------------------------------------------------------------

Usage() {
	echo "Usage: `basename $0` </path/to/folder/>"
	echo
	echo "Example: `basename $0` /biotic/3Tdata/Forbow/085_D/"
	echo
	exit 1
}

SCRIPT=$(python -c "import os; print os.path.abspath('$0')")
SCRIPTSDIR=$(dirname $SCRIPT)
WDIR=$(pwd)

[ "$#" -lt 1 ] && Usage

NCORES=8  ##$(lscpu | egrep 'Core\(s\)' | awk '{print $4}')
DO_POST_RENAME_BZIP2="yes"
VERBOSE="yes"

for arg in "$@" ; do
	pdir=$(python -c "import os; print os.path.abspath('$arg')")
	if [ ! -d "$pdir" ]; then
		echo "* WARNING: could not processs specified Pfile directory: $pdir"
		continue
	fi
	cd $pdir/
	echo "------- running `basename $0` on `pwd`/, starting: `date`"	
	
	
	## skip if no Pnnnnn.7* named files exist
	nPneedRenamed=$(find . -maxdepth 1 -type f -name "P[0-9][0-9][0-9][0-9][0-9].7*" -print | wc -w | bc)
	if [ "$nPneedRenamed" -eq 0 ]; then
		echo " * found no PFILES or param files needing renamed..."
		## find all uncompressed PFILES and compress if nifti file exists and global flag set
		for pfile in $(find . -maxdepth 1 -type f -name "20*_E*_P[0-9][0-9][0-9][0-9][0-9].7" -print); do
			echo " +++ found uncompressed Pfile=$pfile"
			ls -l
			bn=$(basename $pfile)
			if [ -r "${pfile}" -a -r "${bn}.nii.gz" -a "$DO_POST_RENAME_BZIP2" == "yes" ]; then
				#psize=$(du -m $pfile | awk '{x=$1; printf("%d",500+x)}')
				echo " * running cmd: pbzip2 -vf -r -p${NCORES} ${pfile}"
				pbzip2 -vf -r -p${NCORES} ${pfile}
			fi
		done
		continue
	fi
	
	
	## find all BZIP2 compressed Pfiles and create _hdr.txt if none exists already.
	## only decompress PFILES if necessary to run rdgehdr
	for pfile in $(find . -maxdepth 1 -type f -name "P[0-9][0-9][0-9][0-9][0-9].7.bz2" -print); do
		pname=$(basename $pfile '.bz2')
		dc_pfile="$(dirname $pfile)/${pname}"
		phdr="${dc_pfile}_hdr.txt"
		if [ ! -r "${phdr}" ]; then
			pbzip2 -dvf -r -p${NCORES} ${pfile}
			if [ ! -r "${dc_pfile}" ]; then
				echo -e "\n*** ERROR: unable to decompress pfile = ${pfile}...\n"
				exit 5
			fi
			echo " ++ running 'rdgehdr' to create header textfile = $phdr"
			rdgehdr ${dc_pfile} >${phdr}
			if [ ! -r "${phdr}" ]; then
				echo -e "\n*** ERROR: unable to create `pwd`/${phdr}, something is wrong...\n"
				exit 7
			fi
			#psize=$(du -m ${dc_pfile} | awk '{x=$1; printf("%d",500+x)}')
			pbzip2 -vf -r -p${NCORES} ${dc_pfile}
		fi
	done
	
	## assume all Pnnnnn.7 named files now have an associated Pnnnnn.7_hdr...
	for phdr in $(find . -maxdepth 1 -type f -name "P[0-9][0-9][0-9][0-9][0-9].7_hdr.txt" -print); do
		pname=$(basename $phdr '_hdr.txt')
		pfile=$(dirname $phdr)/${pname}
		examNumStr=$(cat -v ${phdr} | grep 'Exam number' | head -1 | awk '{x=$3; printf("E%05d",x);}')
		examTS=$(cat -v ${phdr} | grep 'Exam date' | awk '{print $4}' | egrep -o "[0-9]+")
		examDateStr=$(date -d @$examTS +%Y%m%d-%H%M)
		newName="${examDateStr}_${examNumStr}_${pname}"
		echo " ++ renaming files:  $pname --> $newName"
		rename.ul $pname $newName ${pfile}*
		if [ -r "${newName}" -a -r "${newname}.nii.gz" -a "$DO_POST_RENAME_BZIP2" == "yes" ]; then
			#psize=$(du -m $(dirname $phdr)/${newName} | awk '{x=$1; printf("%d",500+x)}')
			pbzip2 -vf -r -p${NCORES} $(dirname $phdr)/${newName}
		fi
	done
	echo "------- completed `basename $0` on `pwd`/, at: `date`"
	echo
	cd $WDIR/
done

exit 0

