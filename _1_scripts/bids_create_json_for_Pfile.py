#!/usr/bin/env python3

import json, os, sys
import argparse
import time
import re
from scipy.io import loadmat	##needed to read muxrecon .mat file header array


VERBOSE=False

#-----------------------------------------------------------------------------------------
def main(json_file, pfile_hdr, pfile_mat, output_json):
	with open(json_file, 'r') as jfile:
		data = json.load(jfile)

	if VERBOSE: 
		print(' -- deleting unnecessary fields...')
	del data["ConversionSoftwareVersion"]
	del data["ConversionSoftware"]
	del data["ImageOrientationPatientDICOM"]
	del data["InPlanePhaseEncodingDirectionDICOM"]
		
	## get specific values from Pfile_hdr
	with open(pfile_hdr,'r') as pfile:
		phdr_data = pfile.readlines()	

	##  SeriesDescription = ...Series description: EPI MUX 3MM RS
	SerDesc = search_list(re.compile(r"^...Series description:"),phdr_data) 
	if SerDesc == "":
		print(' ** ERROR: could not determine Series Description from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE: 
		print(' ++ grepping "Series description"=%s from Pfile_hdr=%s'%(SerDesc,pfile_hdr))
	data["SeriesDescription"] = SerDesc
	
	##  ProtocolName = ...Scan protocol name: Brain-CANDICE DTI NODDI 
	ProtName = search_list(re.compile(r"^...Scan protocol name:"),phdr_data) 
	if ProtName == "":
		print(' ** ERROR: could not determine ProtName from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE: 
		print(' ++ grepping "ProtocolName"=%s from Pfile_hdr=%s'%(ProtName,pfile_hdr))
	data["ProtocolName"] = ProtName.replace(' ','_')
	
	##  ScanningSequence = ...Pulse sequence name: muxarcepi
	ScanSeq = search_list(re.compile(r"^...Pulse sequence name:"),phdr_data) 
	if ScanSeq == "":
		print(' ** ERROR: could not determine ScanningSequence from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "ScanningSequence"=%s from Pfile_hdr=%s'%(ScanSeq,pfile_hdr))
	data["ScanningSequence"] = ScanSeq
	
	##  SequenceVariant = ...PSD name from inside PSD: EPI
	SeqVar = search_list(re.compile(r"^...PSD name from inside PSD:"),phdr_data) 
	if SeqVar == "":
		print(' ** ERROR: could not determine SeqVar from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "SequenceVariant"=%s from Pfile_hdr=%s'%(SeqVar,pfile_hdr))
	data["SequenceVariant"] = SeqVar.replace(' ','_')
	
	##  ScanOptions = ...Imaging options: 67141632 (Multi-Phase, EPI)
	ScanOpts = search_list(re.compile(r"^...Imaging options:"),phdr_data) 
	if ScanOpts == "":
		print(' ** ERROR: could not determine ScanOptions from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "ScanOptions"=%s from Pfile_hdr=%s'%(ScanOpts,pfile_hdr))
	data["ScanOptions"] = ScanOpts.replace('(','\\').replace(')','\\')
	
	##  ImageType = ??
	
	##  SeriesNumber = ...Series number: 16
	SerNum = search_list(re.compile(r"^...Series number:"),phdr_data) 
	if SerNum == "":
		print(' ** ERROR: could not determine SerNum from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "SeriesNumber"=%s from Pfile_hdr=%s'%(SerNum,pfile_hdr))
	data["SeriesNumber"] = int(SerNum)
	
	##  AcquisitionTime = ...Actual series date/time stamp: (1567591340) Wed Sep  4 07:02:20 2019
	AcqTime = search_list(re.compile(r"^...Actual series date/time stamp:"),phdr_data) 
	if AcqTime == "":
		print(' ** ERROR: could not determine AcqTime from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "AcquisitionTime"=%s from Pfile_hdr=%s'%(AcqTime,pfile_hdr))
	tok=AcqTime.split()
	locAcqTime = time.strftime('%H:%M:%S',time.gmtime(int(tok[0].replace('(','').replace(')',''))))
	data["AcquisitionTime"] = locAcqTime
	
	##  AcquisitionNumber = ...Scan acquisition number: 1
	AcqNum = search_list(re.compile(r"^...Scan acquisition number:"),phdr_data) 
	if AcqNum == "":
		print(' ** ERROR: could not determine AcqNum from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "AcquisitionNumber"=%s from Pfile_hdr=%s'%(AcqNum,pfile_hdr))
	data["AcquisitionNumber"] = int(AcqNum)
	
	##  SliceThickness = ...Slice thickness (mm): 3
	SlcThick = search_list(re.compile(r"^...Slice thickness \(mm\):"),phdr_data) 
	if SlcThick == "":
		print(' ** ERROR: could not determine SlcThick from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "SliceThickness"=%s from Pfile_hdr=%s'%(SlcThick,pfile_hdr))
	data["SliceThickness"] = float(SlcThick)
	
	##  SpacingBetweenSlices = ...Spacing between scans (mm): 0
	SpacBtwSlc = search_list(re.compile(r"^...Spacing between scans \(mm\):"),phdr_data) 
	if SpacBtwSlc == "":
		print(' ** ERROR: could not determine SpacBtwSlc from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "SpacingBetweenSlices"=%s from Pfile_hdr=%s'%(SpacBtwSlc,pfile_hdr))
	if int(SpacBtwSlc) == 0:
		SpacBtwSlc = float(data["SliceThickness"])
	data["SpacingBetweenSlices"] = float(SpacBtwSlc)
	
	##  SAR = ...Peak SAR: 1.71983
	# Possible options from Pfile header:  
	# ...Average SAR: 0.139594
	# ...Peak SAR: 1.71983
	# ...Average head SAR: 0.687931
	SAR = search_list(re.compile(r"^...Peak SAR:"),phdr_data) 
	if SAR == "":
		print(' ** ERROR: could not determine SAR from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "SAR"=%s from Pfile_hdr=%s'%(SAR,pfile_hdr))
	data["SAR"] = float(SAR)
	
	##  EchoTime = ...Pulse echo time (usec): 30000
	# Possible matches form Pfile header
	# *...Pulse echo time (usec): 30000  --> used in mux_epi_recon output.mat file.
	# ...TE for first echo: 30000 
	# ...TE for second and later echoes: 100000
	TE = search_list(re.compile(r"^...Pulse echo time \(usec\):"),phdr_data) 
	if TE == "":
		print(' ** ERROR: could not determine TE from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "EchoTime"=%s from Pfile_hdr=%s'%(TE,pfile_hdr))
	data["EchoTime"] = int(TE)/1000.0/1000
	
	##  RepetitionTime = ...Pulse repetition time (usec): 950000
	TR = search_list(re.compile(r"^...Pulse repetition time \(usec\):"),phdr_data) 
	if TR == "":
		print(' ** ERROR: could not determine TR from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "RepetitionTime"=%s from Pfile_hdr=%s'%(TR,pfile_hdr))
	data["RepetitionTime"] = int(TR)/1000.0/1000
	
	##  FlipAngle = ...Flip angle for GRASS scans (deg): 60
	FlipAng = search_list(re.compile(r"^...Flip angle for GRASS scans \(deg\):"),phdr_data) 
	if FlipAng == "":
		print(' ** ERROR: could not determine FlipAngle from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "FlipAngle"=%s from Pfile_hdr=%s'%(FlipAng,pfile_hdr))
	data["FlipAngle"] = int(FlipAng)
	
	##  PhaseEncodingPolarityGE = ...rhdacqctrl bit mask: 8
	PEdir = search_list(re.compile(r"^...rhdacqctrl bit mask:"),phdr_data) 
	if PEdir == "":
		print(' ** ERROR: could not determine ScanSeq from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "PhaseEncodePolarityGE"=%s from Pfile_hdr=%s'%(PEdir,pfile_hdr))
	if int(PEdir) == 8:
		data["PhaseEncodingPolarityGE"] = "Unflipped"
	elif int(PEdir) == 12:
		data["PhaseEncodingPolarityGE"] = "Flipped"
	else:
		data["PhaseEncodingPolarityGE"] = ""
	
	##  CoilString = ...Coil name: 32Ch Head
	CoilStr = search_list(re.compile(r"^...Coil name:"),phdr_data) 
	if CoilStr == "":
		print(' ** ERROR: could not determine CoilString from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "CoilString"=%s from Pfile_hdr=%s'%(CoilStr,pfile_hdr))
	data["CoilString"] = CoilStr
	
	##  PercentPhaseFOV = ...Display field of view - x (mm): 216
	PhaseFOV = search_list(re.compile(r"^...Display field of view - x \(mm\):"),phdr_data) 
	if PhaseFOV == "":
		print(' ** ERROR: could not determine PercentPhaseFOV from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "PercentPhaseFOV"=%s from Pfile_hdr=%s'%(PhaseFOV,pfile_hdr))
	data["PercentPhaseFOV"] = int(PhaseFOV)
	
	##  AcquisitionMatrixPE = ...Image matrix size - x: 72
	AcqMatrixPE = search_list(re.compile(r"^...Image matrix size - x:"),phdr_data) 
	if AcqMatrixPE == "":
		print(' ** ERROR: could not determine AcquisitionMatrixPE from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "AcquisitionMatrixPE"=%s from Pfile_hdr=%s'%(AcqMatrixPE,pfile_hdr))
	data["AcquisitionMatrixPE"] = int(AcqMatrixPE)
	
	##  ReconMatrixPE = ...Image dimension - x: 72
	RecMatrixPE = search_list(re.compile(r"^...Image dimension - x:"),phdr_data) 
	if RecMatrixPE == "":
		print(' ** ERROR: could not determine ReconMatrixPE from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "ReconMatrixPE"=%s from Pfile_hdr=%s'%(RecMatrixPE,pfile_hdr))
	data["ReconMatrixPE"] = int(RecMatrixPE)
	
	##  ParallelReductionFactorInPlane = ...Number of interleaves: 2
	inplaneR = search_list(re.compile(r"^...Number of interleaves:"),phdr_data) 
	if inplaneR == "":
		print(' ** ERROR: could not determine ParallelReductionFactorInPlane from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "ParallelReductionFactorInPlane"=%s from Pfile_hdr=%s'%(inplaneR,pfile_hdr))
	data["ParallelReductionFactorInPlane"] = int(inplaneR)
	
	
	##  EffectiveEchoSpacing = ...Effective echo spacing for EPI: 232
	##
	## For now, let's assume that EffEchoSpacing reported in Pfile is correct.
	## However, EffEchoSpacing =  (echo-spacing * num-echoes?)
	## - from dcm2niix Github:  #define  kEffectiveEchoSpacingGE  0x0043+(0x102C << 16 ) //SS
	## - 
	EffEchoSpac = search_list(re.compile(r"^...Effective echo spacing for EPI:"),phdr_data) 
	if EffEchoSpac == "":
		print(' ** ERROR: could not determine EffectiveEchoSpacing from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "EffectiveEchoSpacing"=%s from Pfile_hdr=%s'%(EffEchoSpac,pfile_hdr))
	data["EffectiveEchoSpacing"] = int(EffEchoSpac)/1000.0/1000
	
	
	##  "TotalReadoutTime": 0.074992,
	##
	## - https://github.com/rordenlab/dcm2niix/tree/master/GE
	##
	## * Fromula: https://github.com/egarza/tutorials/wiki/Echo-Spacing-and-Total-Readout-Time-Philips-Ingenia-3T 
	##   Total readout time (FSL) = (MatrixSizePhase - 1) * EffectiveEchoSpacing
	##
	## https://bids-specification.readthedocs.io/en/latest/04-modality-specific-files/01-magnetic-resonance-imaging-data.html#diffusion-imaging-data
	##
	data["TotalReadoutTime"] = (float(data["ReconMatrixPE"])-1)*float(data["EffectiveEchoSpacing"])
	
	
	##  "PixelBandwidth" = ...Bandwidth per pixel: 0
	PixBW = search_list(re.compile(r"^...Bandwidth per pixel:"),phdr_data) 
	if PixBW == "":
		print(' ** WARNING: could not determine PixelBandwidth from Pfile_hdr=%s'%(pfile_hdr))
		PixBW = 0.0
	elif VERBOSE:
		print(' ++ grepping "PixelBandwidth"=%s from Pfile_hdr=%s'%(PixBW,pfile_hdr))
	if float(PixBW) > 0.0:
		data["PixelBandwidth"] = float(PixBW)/1000.0/1000
	else:
		## grab bandwidth*2/FreqSpacing: eg. 250kHz * 2.0 / 72.0 = 6944.4444
		VarBW = search_list(re.compile(r"^...Variable bandwidth \(kHz\):"), phdr_data)
		if VarBW == "":
			print(' ** WARNING: could not determine Variable Bandwidth from Pfile_hdr=%s'%(pfile_hdr))
			data["PixelBandwidth"] = 0.0
		else:
			data["PixelBandwidth"] = float(VarBW)*2000/data["ReconMatrixPE"]
			if VERBOSE:
				print(' ++ grepping "Variable Bandwidth"=%s from Pfile_hdr=%s'%(VarBW,pfile_hdr))
	
	## *********************************************************************************
	##  determine Slice-Timing 
	##  - https://cni.stanford.edu/wiki/MUX_EPI
	## *********************************************************************************
	if os.path.exists(pfile_mat):
		mat = loadmat(pfile_mat, squeeze_me=True, struct_as_record=False,variable_names=['p'])
		acq_order = mat['p'].acq_order				# 'interleaved'
		mux_factor = mat['p'].mux					# 3
		sl_order = mat['p'].sl_acq_order.tolist()	# [1, 10, 2, 11, 3, 12, 4, 13, 5, 14, 6, 15, 7, 16, 8, 17, 9]
		num_unmux_sl = mat['p'].num_unmuxed_slices	# 51
		if len(sl_order)*mux_factor != num_unmux_sl:
			print('*** ERROR: num_slices*mux_factor does not match num_unmuxed_slices...')
			sys.exit()
		#sl_start_loc = mat['p'].start_loc			# -9.2566
		data['MultibandAccelerationFactor'] = int(mux_factor)	#"MultibandAccelerationFactor": 4,
		data['ParallelAcquisitionTechnique'] = 'grappa' 		#mux_recon_method: '1Dgrappa'
		if mat['p'].scan_orient == 'axial': 
			data['SliceEncodingDirection'] = 'k'
		else:
			print('*** ERROR: slice-scan orientation is not "axial"...')
			sys.exit()
		nslices = len(sl_order)
		tr = float(data['RepetitionTime'])
		mux_slice_acq_order = sl_order
		mux_slice_acq_time = [float(s)/nslices*tr for s in range(nslices)]
		unmux_slice_acq_order = [nslices*m+s for m in range(mux_factor) for s in mux_slice_acq_order]
		slice_time_list = [slice_time for slice_num,slice_time in zip(unmux_slice_acq_order, mux_slice_acq_time*mux_factor)]
		data['SliceTiming'] = slice_time_list
		if VERBOSE:
			slice_time_dict = {slice_num:slice_time for slice_num,slice_time in zip(unmux_slice_acq_order, mux_slice_acq_time*mux_factor)}
			print(' + SliceTiming:')
			for slice,acq in sorted(slice_time_dict.items()):
				print("  -  slice %02d acquired at time %.3f sec" % (slice+1,acq))		
	
	## SeriesInstanceUID = ...Study entity unique ID: 1.2.840.113619.6.408.18233914799899109266505961543300383295
	SeriesUID = search_list(re.compile(r"^...Study entity unique ID:"),phdr_data) 
	if SeriesUID == "":
		print(' ** ERROR: could not determine SeriesInstanceUID from Pfile_hdr=%s'%(pfile_hdr))
		sys.exit()
	elif VERBOSE:
		print(' ++ grepping "SeriesInstanceUID"=%s from Pfile_hdr=%s'%(SeriesUID,pfile_hdr))
	data["SeriesInstanceUID"] = SeriesUID.lstrip().rstrip().rstrip('\x00')
	### grep -rl '\u0000' * | xargs sed -i'' 's|\\u0000||g'
	
	if VERBOSE:	
		print('\n + json file after modifications:')
		print(json.dumps(data, indent=4))
		print('len(data)=%d'%(len(data)))
	
	## Write out complete .json file
	with open(output_json,'w') as jfile :
		json.dump(data, jfile, indent=4)	
	
	print(' --- saving pfile parameters to file = %s' % output_json)
	return


#-----------------------------------------------------------------------------------------
def search_list(pattern, phdr_list):
	val=""
	for row in phdr_list:
		if pattern.search(row) != None:
			if VERBOSE: print(row.rstrip())
			val=row[row.find(':')+1:].lstrip().rstrip()
			break		
	return val


#-----------------------------------------------------------------------------------------
class ArgumentParser(argparse.ArgumentParser):
    def __init__(self):
        super(ArgumentParser, self).__init__()
        self.description = '''Modifies input PE-matching DWI.json values with appropriate data from specified Pfile_hdr.txt'''
        self.add_argument('-j', dest='json_file', required=True, type=str, metavar='INPUT_JSON_FILE_PATH', help='/path/to/dwi_PE-matching.json')
        self.add_argument('-p', dest='pfile_hdr', required=True, type=str, metavar='PFILE_HDR_FILE_PATH', help='/path/to/Pfile_hdr.txt')
        self.add_argument('-m', dest='pfile_mat', required=True, type=str, metavar='PFILE_MAT_FILE_PATH', help='/path/to/Pfile.mat')
        self.add_argument('-o', dest='output_json', required=True, type=str, metavar='OUTPUT_JSON_PATH', help='/path/to/rs_pe-matching.json')
        self.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
		
		
#-----------------------------------------------------------------------------------------
if __name__ == '__main__':
	args = ArgumentParser().parse_args()
	if args.verbose:
		VERBOSE=True
		print(' +++ Running %s' % os.path.basename(sys.argv[0]))
	main(args.json_file, args.pfile_hdr, args.pfile_mat, args.output_json)
	sys.exit(0)

## ${SCRIPTSDIR}/bids_create_json_for_Pfile.py -j $tmpJson -p $pfileHDR -m $pfileMAT -o $outJsonFile


###########################################################################################
### Pfile.mat from Matlab
###########################################################################################
# p = 
# 
#   struct with fields:
# 
#                        pfile: '20190904-0612_E14874_P17408.7'
#                       FE_DIM: 1
#                       PE_DIM: 2
#                       EC_DIM: 3
#                       SL_DIM: 4
#                        C_DIM: 5
#                        T_DIM: 6
#                       KZ_DIM: 7
#                 KEEP_ORIG_SZ: 1
#                     TOP_DOWN: 0
#                   CENTER_OUT: 1
#                    BOTTOM_UP: 2
#          NUM_COE_PER_KY_LINE: 2
#            ONLY_ACTUAL_COILS: 0
#     ACTUAL_AND_VIRTUAL_COILS: 1
#           ONLY_VIRTUAL_COILS: 2
#                 pfile_header: [1x1 struct]
#                   num_echoes: 1
#                    num_coils: 32
#                   num_passes: 506
#                    inplane_R: 2
#                   partial_ky: 0
#                    dacq_ctrl: 8
#                kissoff_views: 1
#                         vrgf: 1
#                      nx_pres: 72
#                        kydir: 2
#                      ny_pres: 72
#                          fov: 216
#                      slthick: 3
#                    slspacing: 0
#                    start_loc: -9.2566
#                   num_slices: 17
#                  mux_excited: 3
#           num_unmuxed_slices: 51
#                  use_gzblips: 1
#                num_mux_cycle: 2
#                       sldist: 51
#              cap_random_flag: 0
#                       swappf: 1
#                    transpose: 0
#                     rotation: 3
#                   flip_angle: 60
#              slice_shuffling: 0
#                     extra_tr: 0
#                         pseq: 0
#                           r1: 13
#                           r2: 15
#                           te: 30000
#                      epiecho: 1
#           dynamic_phase_corr: 0
#                 fermi_radius: 36
#                  fermi_width: 10
#                    fermi_ecc: 1
#              fermi_rad_scale: 1
#            fermi_width_scale: 1
#               rds_data_order: 0
#                    acq_order: 'interleaved'
#                 sl_acq_order: [17x1 double]
#                 internal_cal: 1
#                  descend_acq: 0
#                  scan_orient: 'axial'
#                multiband_TBW: 1
#                 multiband_pw: 4
#                    isdifscan: 0
#                        next2: 1
#                    dftz_type: 'inverse'
#            cap_fov_shift_cal: -3
#                  mux_encoded: 3
#      recon_cal_use_mux_extra: 0
#                          mux: 3
#                  cal_gzblips: 0
#                        caipi: 1
#                      mica_br: 0
#                    mica_rand: 0
#         mica_perturbed_caipi: 0
#                 mica_poisson: 0
#               cap_blip_start: 2
#                 cap_blip_inc: 2
#                cap_fov_shift: -4
#                       md_ecc: 0
#          ref_for_default_ecc: 'ref.h5'
#                skip_ref_pass: 0
#                  cap_get_ecc: 0
#               coil_noise_std: [1x32 double]
#                     edr_flag: 0
#                        ychop: 1
#            decode_each_slice: 0
#                coil_compress: 0
#                       frames: []
#                        coils: []
#                       echoes: []
#              cal_dat_tpoints: [4 5 6]
#                   zpad_image: 0
#             mux_recon_method: '1Dgrappa'
#                return_sos_im: 1
#              sense_coil_comb: 0
#                grappa_domain: 'imSpace'
#                inplane_kersz: [1x1 struct]
#                inplane_acssz: [1x1 struct]
#          mux_kersz_1d_grappa: [1x1 struct]
#          mux_acssz_1d_grappa: [1x1 struct]
#                zpad_1Dgrappa: 'normalZpad'
#              output_ecc_data: 0
#                      add_vcc: 0
#                        debug: 1
#                  apply_fermi: 0
#                 notch_thresh: 0
#              slices_to_recon: [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17]
#            nt_to_recon_total: 506
#                  nt_to_recon: 500
#              tpoints_to_load: [1x506 double]
#              first_scan_pass: 7
#                     cal_flag: 1
#                      RT_flag: 0
#                      use_GPU: 0
#                    precision: 'double'
#                 calc_snr_pmr: 0
#              calc_sense_rsnr: 0
#                  whiten_type: 'coil_noise_std'
#                 calc_snr_amr: 0
###########################################################################################
