#!/bin/bash

Usage() {
	echo
	echo "Usage: `basename $0` <ssid_sess>"
	echo
	echo "Example: `basename $0` 012_C "
	echo
	exit 1
}


#--------------------------------------------------------------------------------------------------------------
create_mux2nii_script() {
	if [ "$#" -ne 1 -o ! -d "$1" ]; then
		echo "*** ERROR: incorrect usage of function ${0}, exiting..."
		exit 2
	fi
	cat -v >${1}/mux2nii_PBIL.m <<EOF
function mux2nii_PBIL(pfiledir)
% process all p-files within specified "pfiledir"
addpath('$SWDIR/3T_Matlab_Tools_2018/3T_tools');
addpath('$SWDIR/3T_Matlab_Tools_2018/NIfTI');
addpath(genpath('$SWDIR/3T_Matlab_Tools_2018/orchestra-${PFILE_VERS}'));
addpath(genpath('$SWDIR/3T_Matlab_Tools_2018/mux_epi_recon-${PFILE_VERS}'));
if (exist(getenv('FSLDIR'), 'dir') ~= 7),
	fprintf(2,'\n*** ERROR: cannot locate environment FSLDIR in mux2nii_PBIL(),...exiting.\n\n');
	exit(1);
end
if(~exist('pfiledir', 'var')),  %% Locate all p-files within specified "pfiledir"
	pfiledir = pwd;
end
pfiles = dir(strcat(pfiledir,'/201*_E*_P*.7'));	%pattern: "20180509-1336_E09625_P44544.7"
numfiles = length(pfiles);
if(numfiles == 0), 
	fprintf(2,'\n*** ERROR: failed to find any pfiles in specified directory: %s\n\n',pfiledir); 
	exit(1);
end
for ind = 1:numfiles,			%% individually process each pfile
	pfile = pfiles(ind).name;
	if exist([pfile,'.nii'], 'file') == 2 || exist([pfile,'.nii.gz'], 'file') == 2,
		fprintf(1,'\n * WARNING: found existing NIFTI file for pfile=[%s], skipping...\n\n',pfile);
		fprintf(2,'\n * WARNING: found existing NIFTI file for pfile=[%s], skipping...\n\n',pfile);
		continue
	end
	disp(['-- converting pfile: ',pfile]);
	try,
		mux_epi_main_offline(pfile,strcat(pfile,'.mat'));
	catch,
		fprintf(1,'\n*** ERROR: in attempting to run mux_epi_main_offline(pfile,pfile.mat), skipping pfile=[%s]\n\n',pfile);
		fprintf(2,'\n*** ERROR: in attempting to run mux_epi_main_offline(pfile,pfile.mat), skipping pfile=[%s]\n\n',pfile);
		continue
    end
	mux2nii(strcat(pfile,'.mat'));
end
exit(0);

EOF

}
#--------------------------------------------------------------------------------------------------------------


SCRIPT=$(python -c "from os.path import abspath; print(abspath('$0'))")
SCRIPTSDIR=$(dirname $SCRIPT)
if [[ -z "${RAW_DATA_DIR}" ]]; then
	source "$SCRIPTSDIR/FORBOW_SetupEnvironment.sh"
fi

FORCE_OVERWRITE="no"
if [ "$1" == "-f" ]; then
	FORCE_OVERWRITE="yes"
	shift ;
	echo " ** NOTE: Enabling FORCE_OVERWRITE mode..."
fi


############ Setup mux2nii conversion libraries Pfiles ##################################
PFILE_VERS="DV26" ## assume DV26 for all mux-fMRI_Pfiles collected after 2018/01/01

D=$(date +%Y%m%d-%H%M)
FORCE_DECOMPRESS_PRE_CONVERSION="yes"
FORCE_COMPRESS_POST_CONVERSION="yes"
DO_PFILE_CONVERSION="yes"
DO_RPI_REORIENT="yes"
DO_LRFLIP="yes"
DO_NIFTI_COMPRESSION="yes"
PFILE_PATTERN="20??????-????_E*_P*.7"
MATLAB_exe=$(which matlab)
NCORES=${NUM_PHYSICAL_CORES}

[ "$#" -lt 1 ] && Usage

SUBJECTS="$@"
for Subj in $SUBJECTS ; do 
	
	S=$(basename $Subj)
	SDIR=$(find ${RAW_DATA_DIR} -maxdepth 1 -type d -name "${S}_20??????" -print)
	if [ ! -d "$SDIR" ]; then
		echo "*** ERROR: cannot find rawdata directory for $S in ${RAW_DATA_DIR}/"
		continue
	fi
	PFILEDIR="${SDIR}/RS"
	if [ ! -d "$PFILEDIR" ]; then
                echo "*** ERROR: cannot find PFILE directory = $PFILEDIR/"
                continue
        fi
	#mkdir -p ${PFILEDIR}/
	#if [ "`/bin/ls ${PFILEDIR}/* | wc -w | bc`" -lt 12 ]; then
	#	if [ "`ssh biotic@lauterbur /bin/ls -f /biotic/data/Forbow/$S/* | wc -w`" -ge 12 ]; then
	#		echo "+++ transferring PFILES from biotic@lauterbur:/biotic/data/Forbow/${S}/"
	#		rsync -rptl biotic@lauterbur:/biotic/data/Forbow/${S}/* ${PFILEDIR}/
	#	else
	#		echo "*** ERROR: could not locate enough PFILES in folder = $PFILEDIR/"
	#		continue
	#	fi
	#fi
	cd ${PFILEDIR}/
	echo " ----- running `basename $0` on `pwd`/, starting `date`"
	
	
	nonMuxFolder="NON_MUX_PFILES"
	
	## uncompress pfiles, but only if each unzipped Pfile DOES NOT currently exist
	## (to prevent overwriting any expensive data...)
	if [ "$FORCE_DECOMPRESS_PRE_CONVERSION" == "yes" ] ; then
		echo " ++ decompressing PFILES"
		for p in $(find . -type f -name "${PFILE_PATTERN}.gz" -print); do
			pfile=$(basename $p .gz)
			if [[ -r "${pfile}" ]]; then
				echo "* WARNING: found both gzipped and unzipped Pfile = $pfile, skipping decompression and mux2nii conversion..."
				continue
			fi
			if [ "`${FSLDIR}/bin/imtest ${pfile}.nii`" -eq 1 ]; then
				if [ "$FORCE_OVERWRITE" == "no" ]; then
					echo " * Warning: nifti file already exists for Pfile=$pfile, FORCE_OVERWRITE is disabled, skipping decompression and mux2nii conversion..."
					continue
				fi
				echo " * FORCE_OVERWRITE enabled, deleting ${pfile}.nii"
				${FSLDIR}/bin/imrm ${pfile}.nii
			fi
			if [[ ! -r "${pfile}_param.dat" ]] || [[ ! -r "${pfile}_ref.h5" ]] || [[ ! -r "${pfile}_vrgf.dat" ]]; then
				echo "* Note: could not find all required .dat files needed to convert Pfile=$pfile, skipping decompression and mux2nii conversion..."
				continue
			fi	
			pigz -v -p ${NCORES} $p
		done
		for p in $(find . -type f -name "${PFILE_PATTERN}.bz2" -print); do
			pfile=$(basename $p .bz2)
			if [[ -r "${pfile}" ]]; then
				echo "* WARNING: found both bzipped and unzipped Pfile = $pfile, skipping decompression and mux2nii conversion..."
				continue
			fi
			if [ "`${FSLDIR}/bin/imtest ${pfile}.nii`" -eq 1 ]; then
				if [ "$FORCE_OVERWRITE" == "no" ]; then
					echo " * Warning: nifti file already exists for Pfile=$pfile, FORCE_OVERWRITE is disabled, skipping decompression and mux2nii conversion..."
					continue
				fi
				echo " * FORCE_OVERWRITE enabled, deleting ${pfile}.nii"
				${FSLDIR}/bin/imrm ${pfile}.nii
			fi
			if [[ ! -r "${pfile}_param.dat" ]] || [[ ! -r "${pfile}_ref.h5" ]] || [[ ! -r "${pfile}_vrgf.dat" ]]; then
				echo "* Note: could not find all required .dat files needed to convert Pfile=$pfile, skipping decompression and mux2nii conversion..."
				continue
			fi	
			pbzip2 -dv -p${NCORES} $p
		done
	fi
	
	if [ "$DO_PFILE_CONVERSION" == "yes" ]; then
		pfiles=$(find . -type f -name "${PFILE_PATTERN}" -print)
		numPfiles2Convert=$(echo "${pfiles}" | wc -w)
		if [ "${numPfiles2Convert}" -gt 0 ]; then
			echo " ++ beginning matlab 'mux2nii_PBIL()' on the following P-files, starting `date`:"
			/bin/ls -l ${pfiles}
			echo
			ST=${SECONDS}
			create_mux2nii_script $(pwd)
			margs="-nosplash -nojvm -nodesktop -nodisplay"
			${MATLAB_exe} ${margs} -r "try, mux2nii_PBIL('`pwd`'); catch, exit(1); end; exit(0);" </dev/null >mux2nii-${D}.stdout 2>mux2nii-${D}.stderr
			ET=$(($SECONDS - $ST))
			if [ "$?" -ne "0" ]; then
				echo "*** ERROR: mux2nii_PBIL() did not complete correctly, exit code=$?"
			fi
			echo " ++ finished mux2nii_PBIL() on `date`, elapsed time: $(($ET/60)) min, $(($ET%60)) sec"
		else
			echo
			echo " * WARNING: could not find any P-files to convert..."
			echo
		fi
	fi
	
	## always compress NIFTI files - overwrite any pre-existing .nii.gz
	for p in $(find . -type f -name "${PFILE_PATTERN}.nii" -print); do
		pigz -vf -p ${NCORES} $p
	done
	
	## always create timestamp file for each Pfile.
	echo " ++ creating datetimestamp files"
	for nii in $(imglob "${PFILE_PATTERN}.nii.gz"); do
		p=$(basename ${nii} .nii.gz)
		phdr="${p}_hdr.txt"
		tsFile="${p}_timestamp.txt"
		if [[ -r "${phdr}" ]]; then				#already have hdr info, don't need Pfile
			ts=$(cat -v ${phdr} | grep -e '...Actual series date/time' | awk '{print $5}' | egrep -o '[0-9]+' | bc)
			pfile_ts=$(date -r $ts +%Y%m%d%H%M%S)
			echo "$pfile_ts" >${p}_timestamp.txt
		elif [ ! -r "$tsFile" -a -r "$p" ]; then	#don't have a hdr file, need pfile to create timestamp file
			pfile_ts=$(date -r `stat -f "%m" $p` +%Y%m%d%H%M%S)
			echo "$pfile_ts" >${p}_timestamp.txt
		fi
		if [ ! -r "$tsFile" ]; then 
			echo " * Warning: cannot find Pfile=$p to create timestamp file..."
			continue
		fi
		echo " ++ found timestamp file: $p --> $nii --> $pfile_ts"
	done
	
	## re-orient NIFTI files to RPI
	if [ "$DO_RPI_REORIENT" == "yes" ]; then
		for p in $(find . -type f -name "${PFILE_PATTERN}.nii.gz" -print); do
			o=$(3dinfo -orient ${p} | awk '{print $1}')
			if [ "$o" != "RPI" ] ; then
				if [ "$DO_LRFLIP" == "yes" ] ; then
					np="$(basename $p .nii.gz)_origLR.nii.gz"
					${FSLDIR}/bin/immv $p $np
					3dLRflip -X -prefix $p $np
					${FSLDIR}/bin/imrm $np
					echo " ++ flipping LR for file = $p"
				fi
				np="$(basename $p .nii.gz)_${o}.nii.gz"
				${FSLDIR}/bin/immv $p $np
				3dresample -orient RPI -input $np -prefix $p
				${FSLDIR}/bin/imrm $np
				echo " ++ reorienting $p to RPI"
			fi
		done
	fi
	
	## compress Pfiles if requested, there were no conversion errors, Pfile NIFTIs exist, and not already compressed.
	if [ "$FORCE_COMPRESS_POST_CONVERSION" == "yes" ]; then
		if [ x"`cat -v mux2nii-${D}.stderr`" != "x" ]; then
			echo " * Warning: possible errors from mux2nii conversion, skipping compression after conversion to allow for error checking..."
		else
			for p in $(find . -type f -name "${PFILE_PATTERN}" -print); do
				if [ "`${FSLDIR}/bin/imtest ${p}.nii`" -eq 0 ]; then
					echo " * WARNING: could not find NIFTI file for Pfile=$p, skipping compression"
				elif [ -r "${p}.bz2" ]; then
					echo "*** ERROR: found existing compressed version of $p, skipping bzip2 compression..."
				else
					pbzip2 -v -p${NCORES} ${p}
				fi
			done
		fi
	fi
	
	echo " ----- completed running `basename $0` on `pwd`/, at `date`"
	
done

exit 0

