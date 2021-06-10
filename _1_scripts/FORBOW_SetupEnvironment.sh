#!/bin/bash


##### MAIN PATH to MAIN PROJECT ########

# MAIN PATH to CCOMS_2016 data ###########
if [ "`uname -n`" == "Aoraki.local" ]; then
	PROJECT_DIR="$HOME/work/uher/FORBOW"
	export MATLABDIR="/Applications/MATLAB_R2018b.app"
	export ANTSPATH=$HOME/apps/ANTs/bin
	export PATH="$ANTSPATH:$PATH"
else
	PROJECT_DIR="/shared/uher/FORBOW"
	if [ "`uname -s`" == "Darwin" ]; then
		export MATLABDIR="/Applications/MATLAB_R2016b.app"
	else
		export MATLABDIR="/usr/local/MATLAB/R2016b"
	fi
	export ANTSPATH=/opt/ANTs/bin
	export PATH="$ANTSPATH:$PATH"
fi

RAW_DATA_DIR="${PROJECT_DIR}/rawdata"
ANALYSIS_DIR="${PROJECT_DIR}/NIHPD_analysis"
SWDIR="${RAW_DATA_DIR}/_0_software"
export PATH="${SWDIR}/bin:${PATH}"
MNI152DIR="${RAW_DATA_DIR}/_0_software/MNI152"


##### Setup FSL ####################
if [[ -z "$FSLDIR" ]]; then
	if [ "$FSLDIR" != "/usr/local/fsl-6.0.1" ]; then
		export FSLDIR="/usr/local/fsl-6.0.1"
		source ${FSLDIR}/etc/fslconf/fsl.sh
		export PATH="${FSLDIR}/bin:${PATH}" 
		export FSLOUTPUTTYPE="NIFTI_GZ"
	fi
fi
##### Setup AFNI ####################
if [[ -z "$AFNI_DIR" ]]; then
	export AFNI_DIR="/usr/local/afni"
	export PATH="${AFNI_DIR}:${PATH}"
fi
##### Setup MATLAB ####################
if [[ -z "$MATLABDIR" ]]; then
	if [ "`uname -s`" == "Darwin" ]; then
		export MATLABDIR="/Applications/MATLAB_R2016b.app"
	elif [ "`uname -s`" == "Linux" ]; then
		export MATLABDIR="/usr/local/MATLAB/R2016b"
	fi
	export PATH="${MATLABDIR}/bin:${PATH}"
	if [[ -z "$MATLABPATH" ]]; then
		export MATLABPATH="${FSLDIR}/etc/matlab"
	else
		export MATLABPATH="${FSLDIR}/etc/matlab:${MATLABPATH}"
	fi
fi

#### Set OMP variables for parallel aware tools [eddy_openmp, etc] ##############
NUM_PHYSICAL_CORES="4"
if [ "`uname -s`" == "Darwin" ]; then
	NUM_PHYSICAL_CORES=2 ##$(sysctl -n hw.physicalcpu)
elif [ "`uname -s`" == "Linux" ]; then 
	NUM_PHYSICAL_CORES=4 ##$(grep -m 1 'cpu cores' /proc/cpuinfo | awk {'print $4'})
fi
export OMP_NUM_THREADS="$NUM_PHYSICAL_CORES"	


create_1x3_overlay(){
	echo " ++ creating 1x3 overlay from ${1} -> ${2} using lut=${3}"
	slicer ${1} -l ${3} -s 2 -a ${2}
}
create_1x3_RedEdge_overlay(){
	echo " ++ creating 1x3 red-edged overlay from ${1} on ${2} -> ${3}"
	slicer ${1} ${2} -s 2 -a ${3}
}
create_9x1_slicer(){
	echo " ++ creating 9x1 SingleImage from $1 -> $2"
	outImg="$2"
	if [ "${outImg:${#outImg}-4:4}" != ".png" ]; then outImg="${outImg}.png"; fi
	tmp=$(${FSLDIR}/bin/tmpnam)
	slicer ${1} -s 2 -u -x 0.4 ${tmp}_x1.png -x 0.5 ${tmp}_x2.png -x 0.6 ${tmp}_x3.png -y 0.4 ${tmp}_y1.png -y 0.5 ${tmp}_y2.png -y 0.6 ${tmp}_y3.png -z 0.4 ${tmp}_z1.png -z 0.5 ${tmp}_z2.png -z 0.6 ${tmp}_z3.png 
	opts="${tmp}_x1.png + ${tmp}_x2.png + ${tmp}_x3.png + ${tmp}_y1.png + ${tmp}_y2.png + ${tmp}_y3.png + ${tmp}_z1.png + ${tmp}_z2.png + ${tmp}_z3.png"
	pngappend ${opts} ${outImg}
	rm -f ${tmp}*.png
}
create_9x1_midbrain_overlay(){
	echo " ++ creating 9x1 overlay from $1 -> $2 using lut=$3"
	tmp=$(${FSLDIR}/bin/tmpnam)
	slicer ${1} -l ${3} -s 2 -u -x 0.4 ${tmp}_x1.png -x 0.5 ${tmp}_x2.png -x 0.6 ${tmp}_x3.png -y 0.45 ${tmp}_y1.png -y 0.5 ${tmp}_y2.png -y 0.55 ${tmp}_y3.png -z 0.35 ${tmp}_z1.png -z 0.4 ${tmp}_z2.png -z 0.45 ${tmp}_z3.png 
	opts="${tmp}_x1.png + ${tmp}_x2.png + ${tmp}_x3.png + ${tmp}_y1.png + ${tmp}_y2.png + ${tmp}_y3.png + ${tmp}_z1.png + ${tmp}_z2.png + ${tmp}_z3.png"
	pngappend ${opts} ${2}
	rm -f ${tmp}*.png
}
create_9x1_midbrain_RedEdge_overlay(){
	echo " ++ creating 9x1 red-edged overlay from $1 on $2 -> $3"
	tmp=$(${FSLDIR}/bin/tmpnam)
	slicer ${1} ${2} -s 2 -u -x 0.4 ${tmp}_x1.png -x 0.5 ${tmp}_x2.png -x 0.6 ${tmp}_x3.png -y 0.45 ${tmp}_y1.png -y 0.5 ${tmp}_y2.png -y 0.55 ${tmp}_y3.png -z 0.35 ${tmp}_z1.png -z 0.4 ${tmp}_z2.png -z 0.45 ${tmp}_z3.png 
	opts="${tmp}_x1.png + ${tmp}_x2.png + ${tmp}_x3.png + ${tmp}_y1.png + ${tmp}_y2.png + ${tmp}_y3.png + ${tmp}_z1.png + ${tmp}_z2.png + ${tmp}_z3.png"
	pngappend ${opts} ${3}
	rm -f ${tmp}*.png
}
create_9x1_overlay(){
	echo " ++ creating 9x1 overlay from $1 -> $2 using lut=$3"
	tmp=$(${FSLDIR}/bin/tmpnam)
	slicer ${1} -l ${3} -s 2 -u -x 0.4 ${tmp}_x1.png -x 0.5 ${tmp}_x2.png -x 0.6 ${tmp}_x3.png -y 0.4 ${tmp}_y1.png -y 0.5 ${tmp}_y2.png -y 0.6 ${tmp}_y3.png -z 0.4 ${tmp}_z1.png -z 0.5 ${tmp}_z2.png -z 0.6 ${tmp}_z3.png 
	opts="${tmp}_x1.png + ${tmp}_x2.png + ${tmp}_x3.png + ${tmp}_y1.png + ${tmp}_y2.png + ${tmp}_y3.png + ${tmp}_z1.png + ${tmp}_z2.png + ${tmp}_z3.png"
	pngappend ${opts} ${2}
	rm -f ${tmp}*.png
}
create_9x1_RedEdge_overlay(){
	echo " ++ creating 9x1 red-edged overlay of $1 on $2 -> $3"
	tmp=$(${FSLDIR}/bin/tmpnam)
	slicer ${1} ${2} -s 2 -u -x 0.4 ${tmp}_x1.png -x 0.5 ${tmp}_x2.png -x 0.6 ${tmp}_x3.png -y 0.4 ${tmp}_y1.png -y 0.5 ${tmp}_y2.png -y 0.6 ${tmp}_y3.png -z 0.4 ${tmp}_z1.png -z 0.5 ${tmp}_z2.png -z 0.6 ${tmp}_z3.png 
	opts="${tmp}_x1.png + ${tmp}_x2.png + ${tmp}_x3.png + ${tmp}_y1.png + ${tmp}_y2.png + ${tmp}_y3.png + ${tmp}_z1.png + ${tmp}_z2.png + ${tmp}_z3.png"
	pngappend ${opts} ${3}
	rm -f ${tmp}*.png
}
create_12x1_overlay(){
	echo " ++ creating 12x1 overlay from $1 -> $2 using lut=$3"
	tmp=$(${FSLDIR}/bin/tmpnam)
	slicer ${1} -l ${3} -s 2 -u -z -30 ${tmp}_z01.png -z -40 ${tmp}_z02.png -z -50 ${tmp}_z03.png -z -60 ${tmp}_z04.png -z -70 ${tmp}_z05.png -z -80 ${tmp}_z06.png -z -90 ${tmp}_z07.png -z -100 ${tmp}_z08.png -z -110 ${tmp}_z09.png -z -120 ${tmp}_z10.png -z -130 ${tmp}_z11.png -z -140 ${tmp}_z12.png 
	opts="${tmp}_z01.png + ${tmp}_z02.png + ${tmp}_z03.png + ${tmp}_z04.png + ${tmp}_z05.png + ${tmp}_z06.png + ${tmp}_z07.png + ${tmp}_z08.png + ${tmp}_z09.png + ${tmp}_z10.png + ${tmp}_z11.png + ${tmp}_z12.png"
	pngappend ${opts} ${2}
	rm -f ${tmp}*.png
}
create_12x1_RedEdge_overlay(){
	echo " ++ creating 12x1 red-edged overlay from $1 on $2 -> $3"
	tmp=$(${FSLDIR}/bin/tmpnam)
	slicer ${1} ${2} -s 2 -u -z -30 ${tmp}_z01.png -z -40 ${tmp}_z02.png -z -50 ${tmp}_z03.png -z -60 ${tmp}_z04.png -z -70 ${tmp}_z05.png -z -80 ${tmp}_z06.png -z -90 ${tmp}_z07.png -z -100 ${tmp}_z08.png -z -110 ${tmp}_z09.png -z -120 ${tmp}_z10.png -z -130 ${tmp}_z11.png -z -140 ${tmp}_z12.png 
	opts="${tmp}_z01.png + ${tmp}_z02.png + ${tmp}_z03.png + ${tmp}_z04.png + ${tmp}_z05.png + ${tmp}_z06.png + ${tmp}_z07.png + ${tmp}_z08.png + ${tmp}_z09.png + ${tmp}_z10.png + ${tmp}_z11.png + ${tmp}_z12.png"
	pngappend ${opts} ${3}
	rm -f ${tmp}*.png
}
create_8x3_overlay(){
	echo " ++ creating 8x3 overlay from $1 -> $2 using lut=$3"
	tmp=$(${FSLDIR}/bin/tmpnam)
	for ((i=20; i<136; i+=5)); do 
		n=$(printf '%03d' $i); 
		slicer ${1} -l ${3} -s 2 -u -z -${i} ${tmp}${n}.png; 
	done
	pngappend ${tmp}020.png + ${tmp}025.png + ${tmp}030.png + ${tmp}035.png + ${tmp}040.png + ${tmp}045.png + ${tmp}050.png + ${tmp}055.png - ${tmp}060.png + ${tmp}065.png + ${tmp}070.png + ${tmp}075.png + ${tmp}080.png + ${tmp}085.png + ${tmp}090.png + ${tmp}095.png - ${tmp}100.png + ${tmp}105.png + ${tmp}110.png + ${tmp}115.png + ${tmp}120.png + ${tmp}125.png + ${tmp}130.png + ${tmp}135.png ${2}
	rm -f ${tmp}*.png
}
create_8x3_RedEdge_overlay(){
	echo " ++ creating 8x3 red-edged overlay from $1 on $2 -> $3"
	tmp=$(${FSLDIR}/bin/tmpnam)
	for ((i=20; i<136; i+=5)); do 
		n=$(printf '%03d' $i); 
		slicer ${1} ${2} -s 2 -u -z -${i} ${tmp}${n}.png; 
	done
	pngappend ${tmp}020.png + ${tmp}025.png + ${tmp}030.png + ${tmp}035.png + ${tmp}040.png + ${tmp}045.png + ${tmp}050.png + ${tmp}055.png - ${tmp}060.png + ${tmp}065.png + ${tmp}070.png + ${tmp}075.png + ${tmp}080.png + ${tmp}085.png + ${tmp}090.png + ${tmp}095.png - ${tmp}100.png + ${tmp}105.png + ${tmp}110.png + ${tmp}115.png + ${tmp}120.png + ${tmp}125.png + ${tmp}130.png + ${tmp}135.png ${3}
	rm -f ${tmp}*.png
}

report_eTIV_Scalar(){
awkcmd='BEGIN {i=0;} {if(i==0){a=$1;b=$2;c=$3; i=i+1;} else if(i==1){d=$1;e=$2;f=$3; i=i+1;} else if(i==2){g=$1;h=$2;I=$3; i=i+1;}}
END {det=a*e*I+b*f*g+c*d*h-a*f*h-b*d*I-c*e*g; printf("%f\t%f\n",det,1/det);}'
awk "$awkcmd" $1
}
