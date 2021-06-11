#!/bin/bash
# script for execution of deployed applications
#
# Sets up the MATLAB Runtime environment for the current $ARCH and executes 
# the specified command.
#
Usage(){
	echo
	echo "Usage: `basename $0` <deployedMCRroot\> args"
	echo
	echo "Example: `basename $0 </Applications/MATLAB/MATLAB_Runtime/v91> /path/to/pfile`"
	echo
	exit 1
}
exe_name=$(python -c "import os; print os.path.abspath('$0')")
exe_dir=$(dirname "$0")

echo "------------------------------------------"
[ "x$1" = "x" ] && Usage


echo Setting up environment variables
MCRDIR=$(python -c "import os; print os.path.abspath('$1')") #MCRDIR="$1"
if [ "`uname -s`" == "Darwin" ]; then
	OSTYPE="maci64"
	##export DYLD_LIBRARY_PATH=/Applications/MATLAB/MATLAB_Runtime/v91/runtime/maci64:/Applications/MATLAB/MATLAB_Runtime/v91/sys/os/maci64:/Applications/MATLAB/MATLAB_Runtime/v91/bin/maci64
	LIBVAR="${MCRDIR}/runtime/${OSTYPE}:${MCRDIR}/bin/${OSTYPE}:${MCRDIR}/sys/os/${OSTYPE}:${MCRDIR}/sys/opengl/lib/${OSTYPE}:."
	export DYLD_LIBRARY_PATH="${LIBVAR}:${DYLD_LIBRARY_PATH}"
	echo "-- setting DYLD_LIBRARY_PATH = ${DYLD_LIBRARY_PATH}"
else
	OSTYPE="glnxa64"
	LIBVAR="${MCRDIR}/runtime/${OSTYPE}:${MCRDIR}/bin/${OSTYPE}:${MCRDIR}/sys/os/${OSTYPE}:${MCRDIR}/sys/opengl/lib/${OSTYPE}:."
	export LD_LIBRARY_PATH="${LIBVAR}:${LD_LIBRARY_PATH}"
	echo "-- setting LD_LIBRARY_PATH = ${LD_LIBRARY_PATH}"
fi

shift 1
args=
while [ $# -gt 0 ]; do
  token=$1
  args="${args} \"${token}\"" 
  shift
done

eval "\"${exe_dir}/muxrecon\"" $args


exit 0
