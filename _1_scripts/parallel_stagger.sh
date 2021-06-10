#!/bin/bash
#
#/shared/uher/FORBOW/NIHPD_analysis/_1_scripts/parallel_stagger.sh -j 6 -d 600 -s /shared/uher/FORBOW/NIHPD_analysis/_1_scripts/FSv6_run_pipeline_DWI.sh `/bin/ls -dr 027_? 026_* 025_* 024_* 023_* 022_A 020_A`
#
#
function isInt(){
	val="$1"
	case ${val#[-+]} in
		''|*[!0-9]*) echo 0 ;;
		*) echo 1 ;;
	esac
}

function Usage() {
cat <<EOF

Usage: 
  `basename $0` -j <njobs> -d <delay-secs> -s </path/to/script> [<input1> <input2> ... <inputN>]

Examples :
  `basename $0` -j 4 -d 900 -s /shared/uher/FORBOW/NIHPD_analysis/_1_scripts/FSv6_run_pipeline_DWI.sh 001_A 002_A 003_A 004_A 005_A 006_A

Synopsis: This program will stagger the first NJOBS approximately <delay-secs> apart and run in parallel -j NJOBS.

Compulsory Inputs:
  -s <path/to/script>    absolute path to script you want run in parallel providing each input
  inputs                 inputs to give to scripts, (eg, SSIDs)

Options:
  -d <delay-secs>        amount of seconds to stagger each of the first NJOBS apart from each other (default = 300 secs)
  -j <njobs>             number of jobs to run at once
  -v                     enable verbosity
_______________________________________________________________________________________________
EOF
exit 1
}
SCRIPT=$(python -c "from os.path import abspath; print(abspath('$0'))")
SCRIPTSDIR=$(dirname $SCRIPT)
DSTR=$(date +%Y%m%d-%H%M)
######### DEFAULT PARAMETERS ##################################################
TOTAL_CORES=$(parallel --number-of-cores | awk '{print $1}')
echo "-- total cores found on this computer is $TOTAL_CORES"
DEFAULT_NJOBS=$(echo "0.5*$TOTAL_CORES" | bc)
NJOBS="$DEFAULT_NJOBS"
DELAYSECS=900
PROG2RUN=""
VERBOSE=false
SUBJECTS=""
###############################################################################
# Parse arguments
###############################################################################
while (( $# > 1 )) ; do
	case "$1" in
		"-help")
			Usage
			;;
		"-v")
			VERBOSE=true
			echo " -- enabling VERBOSE"
			shift
			;;
		"-d")
			shift
			DELAYSECS="$1"
			if $VERBOSE ; then echo " -- setting DELAYSECS = $DELAYSECS"; fi
			shift
			;;
		"-j")
			shift
			NJOBS="$1"
			if $VERBOSE ; then echo " -- setting NJOBS = $NJOBS"; fi
			shift
			;;
		"-s")
			shift
			PROG2RUN="$1"
			if $VERBOSE ; then echo " -- setting PROGRAM = $PROG2RUN"; fi
			shift
			;;
		-*)
			echo "*ERROR: Unknown option '$1'"
			exit 2
			break
			;;
		*)
			break
			;;
	esac
done

if [ -z "$NJOBS" ] || [ `isInt $NJOBS` -eq 0 ] || [ "$NJOBS" -lt 1 ] ; then
    echo; echo  "ERROR: NJOBS must be > 1..."; echo
    Usage
fi
if [ "$NJOBS" -gt "$TOTAL_CORES" ] ; then
    echo; echo "ERROR: NJOBS must be less than total cores ($TOTAL_CORES)..."; echo
    Usage
fi
if [ -z "$DELAYSECS" ] || [ `isInt $DELAYSECS` -eq 0 ] || [ "$DELAYSECS" -lt 1 ]; then
    echo; echo  "ERROR: DELAYSECS must be > 1..."; echo
    Usage
fi
if [ -z "$PROG2RUN" ] ; then
    echo; echo  "ERROR: script/program to parallelize not specified..."; echo
    Usage
fi
if [ ! -x "$PROG2RUN" ] ; then
    echo; echo  "ERROR: script/program to parallelize does not have execute permissions..."; echo
    Usage
fi
if $VERBOSE ; then
	echo "-- configuring jobfile to run $NJOBS with $DELAYSECS stagger, starting `date +%Y%m%d-%H%M%S` -------"
fi

if [ "$#" -lt 1 ]; then
	echo "ERROR: did not find any inputs to run with ${PROGNAME}..."
	Usage
fi 
INPUTS="$@"
echo "INPUTS:  $INPUTS "

mkdir -p ${SCRIPTSDIR}/jobs
JOBFILE="$SCRIPTSDIR/jobs/parallel_stagger_commands-${DSTR}.jobs"
rm -f $JOBFILE 2>/dev/null
touch $JOBFILE

i=0
for arg in $INPUTS ; do
	d=""
	if [ "$i" -lt "$NJOBS" ] ; then
		d="sleep $(echo "$i * $DELAYSECS" | bc) ; "
	else
		d="sleep 0 ; "
	fi
	c="$d$PROG2RUN $arg "
	if $VERBOSE ; then echo "$c"; fi
	echo "$c" >>$JOBFILE
	let i=i+1
done

cmd="parallel -j $NJOBS -a $JOBFILE"
if $VERBOSE ; then
	echo
	echo "+++++++++++ JOBFILE +++++++++++++++++++++++++++++++++++++++++"
	cat $JOBFILE
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo 
	echo "parallel command: $cmd"
else
	parallel -j ${NJOBS} -a ${JOBFILE}
	rm -f "${JOBFILE}" 2>/dev/null
fi

exit 0
