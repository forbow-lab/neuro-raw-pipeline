function muxrecon(pfile, outfile, ext_cal, slice, n_vcoils, recon_method, vol_flag, GPU_flag)
%
% Simplified main interface for slice-multiplexing EPI reconstruction. (See mux_epi_main_RT for the full interface.)
%
% Inputs
%   pfile       - Filename of the f.pfile, or the directory containing a pfile.
%   outfile     - The file for saving the results.
%   ext_cal     - External calibration. This can be
%                 (1) Filename of the f.pfile, or the directory containing
%                     a pfile. The corresponding external calibration scan
%                     can be either mux>1 or mux=1.
%                 (2) Calibration data matrix for all muxed slices.
%                     Dim: [FE, PE, Echo, MuxedSlice, Coil, Time(mux phase cycling)]
%                 (3) Empty.
%                 For external calibration, 'ext_cal' must not be empty.
%                 For internal calibration, if 'ext_cal'is empty internal
%                 calibration will be used, otherwise external calibration
%                 will be used.
%   slice       - Index of the (muxed) slice to reconstruct.
%   n_vcoils    - Number of virtual coils for coil compression. No coil
%                 compression if n_vcoils is [] or 0 or the-number-of-physical-coils.
%   recon_method- '1Dgrappa' or numbers other than 1: use 1D-GRAPPA;
%                 'sense' or number 1: use SENSE;
%                 'slice-grappa': use slice-GRAPPA;
%                 'split-slice-grappa': use split-slice-GRAPPA.
%                 Append '_sense1' to use sense coil combination. E.g., '1Dgrappa_sense1'.
%                 Default: 1Dgrappa.
%                 This option will be ignored for MICA-type acquisition (always uses SENSE).
%   save_vols   - True: save each volume of each slice into separate mat file.
%                 If set to false, save all volumes of each slice into one mat file.
%   use_GPU     - True: use GPU if the computer has a compatible graphic card. False: only use CPU

if(~exist('vol_flag', 'var'));    vol_flag = '1';     end
if(~exist('GPU_flag', 'var'));    GPU_flag = '1';     end

fprintf('mux_epi_main_RT("%s", "%s", "%s", %d, [], %d, 0, "%s", 1, 1, 0, %d, %s)\n', pfile, outfile, ext_cal, str2num(slice), ...
    str2num(n_vcoils), recon_method, str2num(vol_flag), GPU_flag);
mux_epi_main_RT(pfile, outfile, ext_cal, str2num(slice), [], str2num(n_vcoils), 0, recon_method, 1, 1, 0, str2num(vol_flag), GPU_flag);


