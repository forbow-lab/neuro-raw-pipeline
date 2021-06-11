function d = mux_epi_main_RT(pfile, outfile, ext_cal, slices, nt_to_recon, n_vcoils, debug, recon_method, apply_fermi, use_homodyne,...
                             notch_thresh, save_vols, use_GPU, pmr_reps, noise_fname)
%
% [d, snr_d] = mux_epi_main_RT(f.pfile, [outfile=[]], [ext_cal=[]], [slices=[]], [nt_to_recon=[]], [n_vcoils=[]], [debug=false], [recon_method='1Dgrappa'],
%                           [apply_fermi=false], [use_homodyne=true], [notch_thresh=0], [save_vols=true], [use_GPU=true], [pmr_reps=0], [noise_fname=[]])
%
% Main program for slice-multiplexing EPI reconstruction.
%
% Inputs
%   pfile       - Filename of the pfile, or the directory containing a pfile.
%   outfile     - The file for saving the results. If outfile is empty, then
%                 the results will not be saved to a file.
%   ext_cal     - External calibration. This can be
%                 (1) Filename of the pfile, or the directory containing
%                     a pfile. The corresponding external calibration scan
%                     can be either mux>1 or mux=1.
%                 (2) Calibration data matrix for all muxed slices.
%                     Dim: [FE, PE, Echo, MuxedSlice, Coil, Time(mux phase cycling)]
%                 (3) Empty.
%                 For external calibration, 'ext_cal' must not be empty.
%                 For internal calibration, if 'ext_cal'is empty internal
%                 calibration will be used, otherwise external calibration
%                 will be used.
%   slices      - Indices of the (muxed) slices to reconstruct. Default: All slices.
%   nt_to_recon - Number of time points to reconstruct, excluding the first
%                 few mux phase cycling time points. Not used if 'pfile' has
%                 sequential acquisition. Default: All time points.
%   n_vcoils    - Number of virtual coils for coil compression. No coil
%                 compression if n_vcoils is [] or 0 or the-number-of-physical-coils.
%   debug       - True: calculate intermediate results and print out messages.
%   recon_method- '1Dgrappa' or numbers other than 1: use 1D-GRAPPA;
%                 'sense' or number 1: use SENSE;
%                 'slice-grappa': use slice-GRAPPA;
%                 'split-slice-grappa': use split-slice-GRAPPA.
%                 Append '_sense1' to use sense coil combination. E.g., '1Dgrappa_sense1'.
%                 Default: 1Dgrappa.
%                 This option will be ignored for MICA-type acquisition (always uses SENSE).
%   apply_fermi - True: apply a circular Fermi filter using the rec.fermi parameters
%                 specified in the p-file header.
%   use_homodyne- True: use homodyne to fill partial-k acquisitions. False will use zero-filling.
%                 If use_homodyne is empty, then the p-file header will be used to determine whether
%                 to use homodyne or zero-filling.
%   notch_thresh- If non-zero, a denotching filter will be applied. Any
%                 point in k-space with a value lower than notch_thresh
%                 will be replaced by an adjacent time point.
%   save_vols   - True: save each volume of each slice into separate mat file.
%                 If set to false, save all volumes of each slice into one mat file.
%   use_GPU     - 1: use all GPU devices if the computer has compatible graphic cards; 0: only use CPU.
%                 Array (e.g. [1,2,3]): list of GPU devices to use for the recon. Default: 1.
%   pmr_reps    - Number of repeats in Pseudo Multiple Replica noise map calculation. 
%                 When not set or set to 0, pmr calculation will not run. Default: 0.
%   noise_fname - Noise scan pfile name.
%
% Outputs
%   d           - Reconstructed images. If p.return_sos_im == true, these are square-
%                 root-of-sum-of-squares (SOS) images; if p.return_sos_im == false,
%                 these are complex images. Dim: [1stImageDimension, 2ndImageDimension,
%                 Slice, Time, Coil(=nc_after_coil_compression, if (recon is GRAPPA-type &&
%                 p.return_sos_im == false && p.sense_coil_comb == false); =1, otherwise)].
%                 1stImageDimension - 2ndImageDimension:
%                             R/L - A/P (Axial and Axial-like Oblique)
%                             R/L - S/I (Coronal and Coronal-like Oblique)
%                             A/P - S/I (Sagittal and Sagittal-like Oblique)

% (c) Kangrong Zhu,     Stanford University     Aug 2012


%% Version info of the code, first 7 characters of the hash of the last git commit
gitcommit = '';
if exist('gitlasthash', 'file') 
    fid = fopen('gitlasthash', 'r');
    gitcommit = textscan(fid, '%7c', 1);
    gitcommit = gitcommit{1};
    fclose(fid);
end

%% Parameters for the reconstruction.

% -- File directories and names.
if ~exist('pfile', 'var') || isempty(pfile)
    error('No pfile specified!');
end

mydir = fileparts(mfilename('fullpath'));
addpath(fullfile(mydir, 'utils'), '-end');
addpath(fullfile(mydir, 'pfile_reading'));
addpath(fullfile(mydir, 'coil_compress'));
addpath(fullfile(mydir, 'disp_and_sim_tools'));
addpath(fullfile(mydir, 'correct_ref'));
addpath(fullfile(mydir, 'correct_vrgf'));
addpath(fullfile(mydir, 'espirit'));
addpath(fullfile(mydir, 'pseudo_multiple_replica'));
addpath(fullfile(mydir, 'poisson_disc_pattern'));
addpath(fullfile(mydir, 'sense_rsnr'));

if(~exist('ext_cal', 'var'));      ext_cal     = [];     end
if(~exist('outfile', 'var'));      outfile     = [];     end
if(~exist('slices', 'var'));       slices      = [];     end
if(~exist('nt_to_recon', 'var'));  nt_to_recon = [];     end
if(~exist('n_vcoils', 'var'));     n_vcoils    = [];     end
if(~exist('debug', 'var'));        debug       = false;  end
if(~exist('recon_method', 'var')); recon_method= [];     end
if(~exist('apply_fermi', 'var'));  apply_fermi = false;  end
if(~exist('use_homodyne', 'var')); use_homodyne = true;  end
if(~exist('notch_thresh', 'var')); notch_thresh = 0;     end
if(~exist('save_vols', 'var'));    save_vols = true;     end
if(~exist('use_GPU', 'var'));      use_GPU   = 1;        end
if(~exist('pmr_reps', 'var'));     pmr_reps  = 0;        end
if(~exist('noise_fname', 'var'));  noise_fname = [];     end

if length(slices) > 1
    error('slices must be a scalar in the mux_epi_main_RT call..\n');
end

if debug && ~isempty(outfile)
    time_used = tic;
end

RT_flag = true;     % Flag for real-time recon

% We do lots of ffts on arrays of the same size, so it's worth letting fftw
% measure the optimal algorithm the first time we run one.
% NOTE: this is broken on Octave versions < 3.6.4.
fftw('planner','measure');
if(exist('octave_config_info', 'builtin') && exist(fullfile(mydir,'fftw_wisdom.txt'), 'file'))
   tpoints = load(fullfile(mydir,'fftw_wisdom.txt'));
   fftw('dwisdom',tpoints.wisdom);
end

[f.pfile, f.ref, f.vrgf, f.refp, f.noise_fname, f.pfile_name] = get_pfile_name(pfile);
if ~isempty(noise_fname) && exist(noise_fname, 'file')
    f.noise_fname = noise_fname;            % set to use the noise scan pfile specified by user
end
    
% -- Save all parameters of the mux epi scan and for the reconstruction in a structure.
if debug
    fprintf('Reading header from pfile %s...\n', f.pfile_name);
end
p = mux_epi_params(f.pfile, slices, nt_to_recon, n_vcoils, debug, recon_method, apply_fermi, use_homodyne, notch_thresh, f.noise_fname, f.refp, ...
                   RT_flag, use_GPU, pmr_reps);

if p.internal_cal && (~isempty(ext_cal))
    if p.debug
        fprintf('Will use external calibration data, although header indicates internal calibration data was acquired.\n');
    end
    p.internal_cal = false;
end

if (~p.internal_cal) && isempty(ext_cal)
    error('ext_cal must not be empty if external calibration is to be used.');
end

% -- External calibration
if ~p.internal_cal && ~isnumeric(ext_cal)                                  % ~isnumeric(ext_cal): Input is pfile name or pfile directory, not calibration data matrix
    [f.ext_cal, f.ref_cal, f.vrgf_cal, f.refp_cal, f.noise_fname_cal] = get_pfile_name(ext_cal);

    % Use external cal ref/vrgf for correcting acc data (to avoid unexpected phase variations between external cal and acc scan)
    f.ref = f.ref_cal;
    f.refp = f.refp_cal;
    f.vrgf = f.vrgf_cal;
end

% -- Check parameters
if ~exist(f.ref, 'file') && (strcmp(p.ref_for_default_ecc, 'ref.h5') || strcmp(p.ref_for_default_ecc, 'ref.dat'))    % The ref file is needed for all EPI scans
    error('Cannot find the ref file!');
end

if ~exist(f.vrgf, 'file') && p.vrgf                                        % The vrgf.dat file is needed only when ramp sampling is on
    error('Cannot find the vrgf.dat file!');
end

if p.md_ecc
    if strcmp(p.ref_for_md_ecc, 'ref pfile') && ~exist(f.refp, 'file') && p.cap_get_ecc
        p.ref_for_md_ecc = 'ecc data';
        if p.debug
            fprintf('\nNo reference scan pfile. Will use ECC data collected by mux sequence in matrix-decoding ghosting correctrion.\n\n');
        end
    end
    if strcmp(p.ref_for_md_ecc, 'ecc data') && ~p.cap_get_ecc && exist(f.refp, 'file')
        p.ref_for_md_ecc = 'ref pfile';
        if p.debug
            fprintf('\nNo ECC data collected by mux sequence. Will use f.reference scan f.pfile in matrix-decoding ghosting correctrion.\n\n');
        end
    end
    if (strcmp(p.ref_for_md_ecc, 'ecc data') && ~p.cap_get_ecc) || (strcmp(p.ref_for_md_ecc, 'ref pfile') && ~exist(f.refp, 'file'))
        p.md_ecc = false;
        if p.debug
            fprintf('\nNo reference scan pfile and no ECC data collected by mux sequence, will not use matrix-decoding ghosting correction.\n\n');
        end
    end
end

if ~exist(f.refp, 'file') && (strcmp(p.ref_for_default_ecc, 'ref pfile') || (p.md_ecc && strcmp(p.ref_for_md_ecc, 'ref pfile')))
    error('Cannot find the reference scan pfile!');
end

if p.debug
    fprintf('Scan Parameters:\n');
    fprintf(' MUX excited %d; ARC %d; Number of mux cycling %d; Slice gap %gmm; Matrix [%d, %d];\n', p.mux_excited, p.inplane_R, p.num_mux_cycle, p.sldist, p.nx_pres, p.ny_pres);
    fprintf(' RF TBW=%g, pw=%gms, flip=%d;\n', p.multiband_TBW, p.multiband_pw, p.flip_angle);
    fprintf(' Acquired %d slices. Slice acquisition order (%s):', p.num_slices, p.acq_order); disp(p.sl_acq_order.');
    fprintf(' Acquired %d time points. ', p.num_passes);
    switch p.acq_order
        case 'interleaved'
            fprintf('Number of time points after the %d ecc and %d mux phase cycling to reconstruct is ....\n', p.cap_get_ecc*p.mux_excited, p.mux_encoded*p.num_mux_cycle);
        case 'sequential'
            fprintf('Number of time points to reconstruct is 1, because each TR acquired one slice and equivalently the data contains only 1 time point.\n');
    end
    fprintf(' Internal calibration?\t %d;\n Gz blip on?\t %d;\n', p.internal_cal, p.use_gzblips);
    fprintf(' MUX encoded %d; MUX factor in recon will be %d;\n', p.mux_encoded, p.mux);
    if p.use_gzblips && (p.mux_excited > 1)
        if p.caipi
            fprintf(' CAIPI acquisition, CAIPI FOV shift for accelerated time points is %d;\n', p.cap_fov_shift);
        end
        if p.mica_br
            fprintf(' Bit-reversed MICA acquisition, temporal seed shift = %d;\n', p.cap_seed_shift);
        end
        if p.mica_rand
            fprintf(' Random MICA acquisition, temporal seed shift = %d;\n', p.cap_seed_shift);
        end
        if p.mica_perturbed_caipi
            fprintf(' Perburbed CAIPI MICA acquisition, temporal seed shift = %d, kz perturb range = %.4f;\n', p.cap_seed_shift, p.cap_kz_rand_pert_range);
        end
        if p.mica_poisson
            fprintf(' Poisson-disc MICA acquisition;\n');
        end
    end
    fprintf(' Ramp sampling on?\t %d;\n Partial ky?\t %d;\n descend_acq?\t %d;\n ychop?\t %d;\n\n', p.vrgf, p.partial_ky, p.descend_acq, p.ychop);

    if(p.rds_data_order)
        fprintf('RDS was used to acquire data; assuming rds data ordering for the raw data\n');
    end
end

% -- External calibration
if ~p.internal_cal && ~isnumeric(ext_cal)
    p_cal = mux_epi_params(f.ext_cal, slices, [], n_vcoils, debug, recon_method, apply_fermi, use_homodyne, notch_thresh, f.noise_fname_cal, f.refp_cal, RT_flag, use_GPU);

    p.ref_for_default_ecc = p_cal.ref_for_default_ecc;
   
    % Sometimes the coil_noise_std in the target scan isn't set correctly.
    % In that case, use the values from the cal scan.
    if (~isfield(p, 'coil_noise_std') || isempty(p.coil_noise_std) || any(isnan(p.coil_noise_std))) && isfield(p_cal, 'coil_noise_std')
        p.coil_noise_std = p_cal.coil_noise_std;
    end

    if (p_cal.mux>1 && p_cal.num_mux_cycle <= 0); error('Ext Cal: num_mux_cycle <= 0.\n'); end
    if (p_cal.mux>1 && p_cal.sldist ~= p.sldist); error('Ext Cal: Band separation mismatch.\n'); end
    if (p_cal.nx_pres ~= p.nx_pres); error('Ext Cal: Prescribed size in FE mismatch.\n'); end
    if (p_cal.ny_pres ~= p.ny_pres); error('Ext Cal: Prescribed size in PE mismatch.\n'); end
    if (p_cal.inplane_R ~= p.inplane_R); disp('Ext Cal: Inplane acceleration factor mismatch.\n'); end
    if (p_cal.mux>1 && p_cal.multiband_TBW ~= p.multiband_TBW); disp('Ext Cal: TBW of RF pulse mismatch.\n'); end
    if (p_cal.mux>1 && p_cal.multiband_pw ~= p.multiband_pw); disp('Ext Cal: Pulse duration of RF pulse mismatch.\n'); end

    switch p_cal.mux_excited
        case 1                                                           % MUX=1 external calibration
            nt_to_recon_cal = [];                                        % Load all time points (including the first p_cal.num_mux_cycle time points)
            p_cal.num_mux_cycle = 1;
        case p.mux_excited                                               % MUX=p.mux_excited external calibration
            nt_to_recon_cal = 0;                                         % Load all mux cycling data (don't load accelerated data)
            if p_cal.num_passes < (p_cal.cap_get_ecc*p_cal.mux_excited + p_cal.mux*p_cal.num_mux_cycle) % For external calibration whose actual number of time points is smaller than the number of time points needed for the specified num_mux_cycle
                p_cal.num_mux_cycle = floor((p_cal.num_passes-p_cal.cap_get_ecc*p_cal.mux_excited)/p_cal.mux); % Keep only fully sampled mux cyclings
            end
        otherwise
            error('Ext Cal: MUX factor must be either 1 or the same as the actual scan.\n');
    end
    p_cal = mux_epi_params_set_tpoints(p_cal, nt_to_recon_cal);

    if strcmp(p_cal.ref_for_default_ecc, 'ref pfile') || (p_cal.md_ecc && strcmp(p_cal.ref_for_md_ecc, 'ref pfile')) || ...
            (ismember(p_cal.mux_recon_method, {'slice-grappa', 'split-slice-grappa'}) && p_cal.slice_grappa_odd_even_fit && strcmp(p_cal.ref_for_md_ecc, 'ref pfile'))
        if ismember(p_cal.mux_recon_method, {'slice-grappa', 'split-slice-grappa'})
            if p_cal.slice_grappa_odd_even_fit
                md_ecc = true;
            else
                md_ecc = false;
            end
        else
            md_ecc = p_cal.md_ecc;
        end
        p_cal.ref_slices = get_ref_slices(f.refp_cal, p_cal.slices_to_recon, p_cal.num_slices, p_cal.mux, md_ecc, p_cal.debug);
    end

    % If external calibration has longer DFT encoding, need to update a few parameter fields
    if (p_cal.mux_excited == p.mux_excited) && (p_cal.mux_encoded > p.mux_encoded)
        p.recon_cal_use_mux_extra = p_cal.recon_cal_use_mux_extra;
        if p.mux ~= p_cal.mux
            p.mux = p_cal.mux;
            p.num_vcoils = p_cal.num_vcoils;
            p.num_unmuxed_slices = p_cal.num_unmuxed_slices;
            if strcmp(p.mux_recon_method, '1Dgrappa')
                p.mux_acssz_1d_grappa.y = p_cal.mux_acssz_1d_grappa.y;
                p.zpad_1Dgrappa = p_cal.zpad_1Dgrappa;
            end
        end
    end
end

% Save calibration type in a variable
if p.internal_cal
    cal_type = 'internal';
else
    if p_cal.mux_excited == 1
        cal_type = 'external_single_band';
    elseif p_cal.mux_excited == p.mux_excited
        cal_type = 'external_multiband';
    end
end

% Slices(Normally ordered indices, not pfile indices) in reference scan corresponding to the slices(Normally ordered indices, not pfile indices) to reconstruct in actual mux scan
if strcmp(p.ref_for_default_ecc, 'ref pfile') || (p.md_ecc && strcmp(p.ref_for_md_ecc, 'ref pfile')) || ...
        (ismember(p.mux_recon_method, {'slice-grappa', 'split-slice-grappa'}) && p.slice_grappa_odd_even_fit && strcmp(p.ref_for_md_ecc, 'ref pfile'))
    if ismember(p.mux_recon_method, {'slice-grappa', 'split-slice-grappa'})
        if p.slice_grappa_odd_even_fit
            md_ecc = true;
        else
            md_ecc = false;
        end
    else
        md_ecc = p.md_ecc;
    end
    p.ref_slices = get_ref_slices(f.refp, p.slices_to_recon, p.num_slices, p.mux, md_ecc, p.debug);
end

%% Call the reconstruction routines.
% Some parameters to deal with pfile reading hit the end of the file while it's still being written
aborted_scan = false;

% We have to actually allocate this below, becuase we don't yet know the timeseries length.
d = [];
pmr = [];
tpoints_step = 1;

for sl = p.slices_to_recon(:)'

    if p.debug
        fprintf('Processing slice %d/%d...\n', find(p.slices_to_recon == sl), length(p.slices_to_recon));
    end

    % -- Slice indices in reference scan
    if strcmp(p.ref_for_default_ecc, 'ref pfile') || (p.md_ecc && strcmp(p.ref_for_md_ecc, 'ref pfile')) || ...
            (ismember(p.mux_recon_method, {'slice-grappa', 'split-slice-grappa'}) && p.slice_grappa_odd_even_fit && strcmp(p.ref_for_md_ecc, 'ref pfile'))
        sl_idx = find(p.slices_to_recon == sl);
        sl_ref = p.ref_slices(sl_idx : length(p.slices_to_recon) : end);
    elseif ~p.internal_cal && (p_cal.mux_excited == 1)                  % if using ref.dat of a single band external calibration scan
        sl_ref = sl + floor((p.mux - 1)/2) * p.num_slices;              % then use the pha_coe of the center band slices for correcting the mux slices
    else
        sl_ref = sl;
    end
 

    % (Re-)Initialize some parameters for each slice
    dat = [];
    d_sl = [];
    p.cal_flag = 1;
    p.num_coils = p.pfile_header.ncoils;   % this parameter is changed in coil compression, which should actually be saved in p.num_vcoils
    p.ccmtx = [];
    p.smap  = [];
    p.solve_mtx = [];
    p.inplane_ker = [];
    if strcmp(p.mux_recon_method, 'slice-grappa') || strcmp(p.mux_recon_method, 'split-slice-grappa')
        p.through_plane_ker = {};
    else
        p.through_plane_ker = [];
    end
    
    if p.internal_cal
        p.cal_dat_tpoints = p.skip_ref_pass + (p.mux*(p.num_mux_cycle-1)+1 : p.mux*p.num_mux_cycle);
        p.tpoints_to_load = p.skip_ref_pass + (1 : (p.mux_excited*p.cap_get_ecc + p.mux_encoded*p.num_mux_cycle));
    else
        p.cal_dat_tpoints = p_cal.skip_ref_pass + (p_cal.mux*(p_cal.num_mux_cycle-1)+1 : p_cal.mux*p_cal.num_mux_cycle);
        p.tpoints_to_load = p_cal.skip_ref_pass + (1 : (p_cal.mux_excited*p_cal.cap_get_ecc + p_cal.mux_encoded*p_cal.num_mux_cycle));
    end

    while p.tpoints_to_load(1) <= p.nt_to_recon_total

        if p.internal_cal || (~p.internal_cal && ~p.cal_flag)

            if p.debug
                fprintf('Loading data from p-file...\n');
            end

            % -- Load the EPI time series data (Already corrected for EPI ghosting and ramp sampling).
            sl_in_pfile = p.sl_acq_order(sl);
            [dat, p, dat_ecc] = epi_load_tseries(f.pfile, f.ref, f.vrgf, f.refp, sl, sl_in_pfile, sl_ref, p); % p_recon will be used as the parameter structure in the reconstruction. Several fields may have been changed in function epi_load_tseries. Parameters many change later if external calibration or coil compression is used. Can't modify values in p because it is reused by multiple slices.

            if isempty(dat)
                aborted_scan = true;
                break;
            end
        end
        
        % -- For debugging: Check the effect of removing the encoding added to one of the simultaneous slices
        if p.debug && p.decode_each_slice
            if (p.mica_br || p.mica_rand) && (p.tpoints_to_load(end) - p.mux_excited*p.cap_get_ecc - p.mux_encoded*p.num_mux_cycle > 0)
                nt_decode = p.tpoints_to_load(end) - p.mux_excited*p.cap_get_ecc - p.mux_encoded*p.num_mux_cycle;
            else
                nt_decode = 1;
            end
            dat_decoded_each_slice = mux_decode_each_slice(dat(:, 1:p.inplane_R:end, :, :, :, p.mux*p.num_mux_cycle+(1:nt_decode)), ...
                p.mux, get_ky_omegaz_us_msk(p, size(dat, p.PE_DIM), nt_decode, true));

            % Check synthesized mux data vs. actual mux data
            if p.use_gzblips
                p.us_msk = get_ky_omegaz_us_msk(p, size(dat, p.PE_DIM), nt_decode, true);
                dat_cal = dat(:, :, 1, 1, :, p.cal_dat_tpoints);
                dat_cal = 1/size(dat_cal, p.T_DIM) * fft(dat_cal, [], p.T_DIM);
                dat_cal_sz = get_dat_sz(dat_cal, p);
                nt_compare = size(dat, p.T_DIM) - p.mux*p.num_mux_cycle;
                ec = 1;
                sl_idx = find(p.slices_to_recon == sl);
                if ~exist('dat_acc', 'var')
                    dat_acc = zeros(dat_cal_sz.x, dat_cal_sz.y, dat_cal_sz.ec, length(p.slices_to_recon), dat_cal_sz.c, nt_compare, 'single');
                end
                dat_acc(:, :, ec, sl_idx, :, :) = dat(:, :, 1, 1, :, p.mux*p.num_mux_cycle+1:end);
                if ~exist('dat_syn', 'var')
                    dat_syn = zeros(dat_cal_sz.x, dat_cal_sz.y, dat_cal_sz.ec, length(p.slices_to_recon), dat_cal_sz.c, nt_compare, 'single');
                end
                for time = 1 : nt_compare
                    if length(p.us_msk) == 1
                        us_msk_this_time = p.us_msk;
                    else
                        us_msk_this_time = p.us_msk(time);
                    end
                    for coil = 1 : dat_cal_sz.c
                        for x = 1 : dat_cal_sz.x
                            tpoints = dat_cal(x, 1:p.inplane_R:end, ec, 1, coil, :);
                            tpoints = reshape(tpoints, [dat_cal_sz.y/p.inplane_R, p.mux]);
                            tpoints = tpoints .* ifftshift(us_msk_this_time.ftz_pha, 2);
                            dat_syn(x, 1:p.inplane_R:end, ec, sl_idx, coil, time) = sum(tpoints, 2);
                        end
                    end
                end
            end
        end

        % -- External calibration
        if ~p.internal_cal && p.cal_flag
            if ~isnumeric(ext_cal)                                           % Input is pfile name or pfile directory
                if p.debug
                    fprintf(' Loading data from external calibration p-file...\n');
                    fprintf('  Number of MUX cycles in external calibration is %d.\n', p_cal.num_mux_cycle);
                end

                switch p_cal.mux_excited
                    case 1                                                   % MUX=1 external calibration
                        sl_cal = sl : p.num_slices : sl+(p.mux-1)*p.num_slices;
                        if p_cal.descend_acq && strcmp(p.dftz_type, 'inverse')% flip the slice order of calibration if single band external cal has descending slices, because caipi encoding in the SMS scan is always ascending
                            sl_cal = flip(sl_cal);
                        end
                        sl_in_pfile_cal = p_cal.sl_acq_order(sl_cal);
                    case p.mux_excited                                       % MUX=p.mux_excited external calibration
                        sl_cal = sl;
                        sl_in_pfile_cal = p_cal.sl_acq_order(sl);
                end
                
                % -- Slice indices in reference scan
                if strcmp(p_cal.ref_for_default_ecc, 'ref pfile') || (p_cal.md_ecc && strcmp(p_cal.ref_for_md_ecc, 'ref pfile')) || ...
                        (ismember(p_cal.mux_recon_method, {'slice-grappa', 'split-slice-grappa'}) && p_cal.slice_grappa_odd_even_fit && strcmp(p_cal.ref_for_md_ecc, 'ref pfile'))
                    switch p_cal.mux_excited
                        case 1
                            sl_ref_cal = repmat(sl + floor((p.mux - 1)/2) * p.num_slices, 1, p.mux); % use the phase correction coef of center band slices for correcting all the slices that will be excited simultaneously in the accelerated timepoints
                        case p.mux_excited
                            sl_idx = find(p_cal.slices_to_recon == sl);
                            sl_ref_cal = p_cal.ref_slices(sl_idx : length(p.slices_to_recon) : end);
                    end
                elseif (p_cal.mux_excited == 1)                                     % if using ref.dat of a single band external calibration scan
                    sl_ref_cal = repmat(sl + floor((p.mux - 1)/2) * p.num_slices, 1, p.mux); % then use the pha_coe of the center band slices for correcting the mux slices    
                else
                    sl_ref_cal = sl;                                                % if using ref.dat of a mux external calibration, then use the pha_coe of the same mux slice from the calibration for phase correcting the accelerated timepoints
                end

                if p_cal.isdifscan
                    p_cal.tpoints_to_load = p_cal.skip_ref_pass + (1 : (p_cal.mux_excited*p_cal.cap_get_ecc + p_cal.mux_encoded*p_cal.num_mux_cycle));
                end
                [dat_cal, p_cal_recon] = epi_load_tseries(f.ext_cal, f.ref_cal, f.vrgf_cal, f.refp_cal, sl_cal, sl_in_pfile_cal, sl_ref_cal, p_cal);

                if isempty(dat_cal)
                    aborted_scan = true;
                    if p.debug
                        warning('External calibration file doesn''t have enough data for calibration.');
                    end
                    break;
                end
            
                if (p_cal_recon.mux_excited == p.mux_excited) && (p_cal_recon.mux_encoded ~= p.mux_encoded)
                    p.cap_fov_shift_cal = p_cal_recon.cap_fov_shift_cal; % p_cal_recon.cap_fov_shift_cal and p_cal_recon.mux_encoded may have been changed in function epi_load_tseries
                    p.mux_encoded = p_cal_recon.mux_encoded;
                end
                
                p_cal_recon_fields = fieldnames(p_cal_recon);
                p_cal_recon_values = struct2cell(p_cal_recon);
                for i = 1:length(p_cal_recon_fields)
                    if ~isfield(p, p_cal_recon_fields{i})
                        p = setfield(p, p_cal_recon_fields{i}, p_cal_recon_values{i});
                    end
                end

                switch p_cal.mux_excited
                    case 1                                                   % MUX=1 external calibration
                        dat_cal_num_mux_cycle = 1;                           % Number of mux cycles in the synthesized multiplexed calibration data
                        if p_cal.inplane_R > 1
                            dat_cal = dat_cal(:, :, :, :, :, p_cal.num_mux_cycle); % Use the last fully sampled time point since all later time points are undersampled.
                        else
                            dat_cal = dat_cal(:, :, :, :, :, p_cal.num_mux_cycle:end); % Exclude the first few non steady state time points
                            dat_cal = mean(dat_cal, p.T_DIM);                % Average all time points for higher SNR
                        end
                        dat_cal = ifftshift(dat_cal, p.SL_DIM);              % Make the middle slice the 0th indexed slice. Dim: [Kx, Ky, Echo, SimultaneousSlice Z(=p.mux), Coil, Time(=1)]
                        dat_cal = mux_dftz(dat_cal, p.SL_DIM, p.cap_fov_shift_cal, p.mux, 'encode');
                        dat_cal = permute(dat_cal, [1,2,3,6,5,4]);           % Dim: [Kx, Ky, Echo, Slice(=1), Coil, Time(=p.mux)]
                    case p.mux_excited                                       % MUX=p.mux_excited external calibration
                        dat_cal_num_mux_cycle = p_cal.num_mux_cycle;
                end
            else                                                             % Input is calibration data matrix
                dat_cal_num_mux_cycle = size(ext_cal, p.T_DIM) / p.mux;
                dat_cal = ext_cal(:, :, :, sl, :, :);
            end

            %dat = dat(:, :, :, :, :, p.mux_encoded*p.num_mux_cycle+1:end);   % Discard internal calibration data, if there were any
            %dat = cat(p.T_DIM, dat_cal, dat);
            dat = dat_cal;
            dat_num_mux_cycle = p.num_mux_cycle;                       % store the original num_mux_cycle of the actual scan
            p.num_mux_cycle = dat_cal_num_mux_cycle;                   % num_mux_cycle must be changed after loading the data from the actual scan and before reconstruction
            p.cal_dat_tpoints = p.mux*(p.num_mux_cycle-1)+1 : p.mux*p.num_mux_cycle; % Time points corresponding to the calibration data
            p.mux_raw_data = dat;                                  % store multiplexed calibration in the parameters for slice-grappa recon
        end

        % -- Coil Compression.
        if p.coil_compress
            if p.debug
                fprintf(' Conducting coil compression...\n');
            end

            [dat, p] = coil_compress(dat, p);
        end

        % -- Recon the data.
        if p.debug
            fprintf(' Reconstructing...\n');
        end
        switch p.mux_recon_method
            case 'sense'
                [dat, p] = mux_epi_process_data_sense(dat, p);
            case '1Dgrappa'
                [dat, p] = mux_epi_process_data_grappa(dat, p);
            case 'slice-grappa'
                [dat, p] = mux_epi_process_data_slice_grappa(dat, p);
            case 'split-slice-grappa'
                [dat, p] = mux_epi_process_data_slice_grappa(dat, p);
            otherwise
                error(['Unknown recon method "' p.mux_recon_method '". Aborting.']);
        end

        % -- Add to results the ECC data collected by the mux sequence
        if p.output_ecc_data && ~p.zpad_image
            if size(dat, p.C_DIM) == 1
                dat_ecc = sos(dat_ecc, p.C_DIM);
            end
            if p.partial_ky
                dat_ecc_sz = get_dat_sz(dat_ecc, p);
                dat_ecc = cat(p.PE_DIM, dat_ecc, zeros(dat_ecc_sz.x, p.ny_pres-dat_ecc_sz.y, dat_ecc_sz.ec, dat_ecc_sz.sl, dat_ecc_sz.c, dat_ecc_sz.t, 'single'));
            end
            dat = cat(p.T_DIM, dat_ecc, dat);
        end

        if p.debug
            fprintf(' Preparing outputs...\n');
        end
        % -- Fix phase-encode direction for pepolar scans
        if bitand(p.dacq_ctrl, 4) % 3rd bit indicates odd-echo phase-flip (i.e., "pepolar")
            dat = flip(dat, p.PE_DIM);
        end

        % -- Fix slice ordering.
        if p.epiecho > 1                                                    % Use the precribed number of slices to find sl_loc
            num_slices = p.num_slices / p.epiecho;
            num_unmuxed_slices = p.num_unmuxed_slices / p.epiecho;
            n_echo = floor((sl-1)/num_slices) + 1;                          % Record the echo number of the reconstructed slice
            sl_orig = sl - (n_echo-1) * num_slices;
        else
            num_slices = p.num_slices;
            num_unmuxed_slices = p.num_unmuxed_slices;
            n_echo = 1;
            sl_orig = sl;
        end            
        if p.descend_acq                                                     % If the multiplexed slices are acquired in descending order. Need this because the unmuxed slices are always in ascending order (i.e., from bottom to top) for the current DFT encoding scheme in the PSD.
            sl_loc = (num_slices*(p.mux-1)+sl_orig) : -num_slices : 1;
        else
            sl_loc = sl_orig : num_slices : num_unmuxed_slices;
        end

        % -- Save the recon results.
        % Note that size(im,6) might be less than the prescribed timeseries length
        % (i.e., p.num_mux_cycle+p.nt_to_recon) if the scan was stopped prematurely.
        if p.zpad_image
            nx_image = p.zpad_size(1);
            ny_image = p.zpad_size(2);
        else
            nx_image = p.nx_pres;
            ny_image = p.ny_pres;
        end
        if p.return_sos_im                                                                      % Return SOS coil-combined magnitude image
            d_size = [nx_image, ny_image, p.num_unmuxed_slices, size(dat, p.T_DIM)];             % Dim: [X, Y, Slice, Time]
        else                                                                                    % Return complex images
            if strcmp(p.mux_recon_method, 'sense')                                              % SENSE recon, return coil-combined complex image
                d_size = [nx_image, ny_image, p.num_unmuxed_slices, size(dat, p.T_DIM)];         % Dim: [X, Y, Slice, Time]
            else                                                                                % GRAPPA recon, return complex single-coil images (when p.sense_coil_comb is false) or coil-combined complex image (when p.sense_coil_comb is true)
                if ~p.sense_coil_comb
                    d_size = [nx_image, ny_image, p.num_unmuxed_slices, size(dat, p.T_DIM), size(dat, p.C_DIM)]; % Dim: [X, Y, Slice, Time, Coil]
                else
                    d_size = [nx_image, ny_image, p.num_unmuxed_slices, size(dat, p.T_DIM)];     % Dim: [X, Y, Slice, Time]
                end
            end
        end
        if p.isdifscan && (p.next2 > 1)
            dat = mean(dat, 3);
        end
        if p.return_sos_im
            dat = reshape(dat, [nx_image, ny_image, p.mux, size(dat, p.T_DIM)]);
        else
            if strcmp(p.mux_recon_method, 'sense')
                dat = reshape(dat, [nx_image, ny_image, p.mux, size(dat, p.T_DIM)]);
            else
                if ~p.sense_coil_comb
                    dat = permute(reshape(dat, [nx_image, ny_image, p.mux, size(dat, p.C_DIM), size(dat, p.T_DIM)]), [1,2,3,5,4]); % Dim: [X, Y, Slice, Time, Coil]
                else
                    dat = reshape(dat, [nx_image, ny_image, p.mux, size(dat, p.T_DIM)]);
                end
            end
        end

        % -- Adjust orientation of reconed images
        if p.transpose
            dat = permute(dat, [2,1,3,4,5]); % Dim: [1stImageDimension, 2ndImageDimension, Slice, Time, Coil(=nc_after_coil_compression if d contains complex images from GRAPPA recon; =1 otherwise)]
        end
        if p.rotation>0
            sz = [size(dat, 1), size(dat, 2), size(dat, 3), size(dat, 4), size(dat, 5)];
            if mod(p.rotation,2)
                d_size([1 2]) = d_size([2 1]);
                sz([1 2]) = sz([2 1]);
            end
            tmp = dat;
            dat = zeros(sz);
            for i=1:sz(3)
                for j=1:sz(4)
                    for k=1:sz(5)
                        dat(:,:,i,j,k) = rot90(tmp(:,:,i,j,k), p.rotation);
                    end
                end
            end
        end

        if p.swappf
            fprintf('Swapping freq/phase dims...\n');
            dat = permute(dat, [2,1,3,4,5]); % Dim: [1stImageDimension, 2ndImageDimension, Slice, Time, Coil(=nc_after_coil_compression if d contains complex images from GRAPPA recon; =1 otherwise)]
            d_size([1 2]) = d_size([2 1]);
        end

        % -- Scale the images according to the receive gains
        % scaling = 10 * 2^((p.r2 - 15) + 0.5 * (p.r1 - 14));
        if p.edr_flag
            scaling = 1/2^16;
        else
            scaling = 100;
        end
        dat = dat * scaling;

        % -- Save results
        if save_vols && (~isempty(dat) && ~isempty(outfile))  % Save 'd' to mat file at each slice, each time point
            d = dat(:, :, :, end, :);
            if exist('octave_config_info', 'builtin')
                save(sprintf([outfile,'_t%04d_s%03d'], p.nt_to_recon, sl-1), 'd', 'sl_loc', 'd_size', '-v7');
            else
                if ~(p.debug && p.decode_each_slice && p.use_gzblips)
                    dat_acc = [];
                    dat_syn = [];
                end
                if p.debug && p.decode_each_slice
                    im_decoded_each_slice = sos(ifft2c(dat_decoded_each_slice), p.C_DIM);
                else
                    im_decoded_each_slice = [];
                end

                if p.internal_cal || (~p.internal_cal && ~p.cal_flag)
                    outfile_vols = sprintf([outfile,'_t%04d_s%03d'], p.nt_to_recon, sl-1);
                else    % if using ext_cal, the first volume (t0000) will be calibration, all accelerated time points start from t0001
                    outfile_vols = sprintf([outfile,'_t%04d_s%03d'], 0, sl-1);
                end
                save(outfile_vols, 'd', 'sl_loc', 'd_size', 'n_echo', 'cal_type', 'ext_cal', 'pfile', 'gitcommit', '-v7');
            end
        elseif ~save_vols && ~isempty(dat)   % Save all time points of each slice to 'd_sl'
                d_sl = cat(4, d_sl, dat);
        end
        
        if p.internal_cal || (~p.internal_cal && ~p.cal_flag)
            if p.tpoints_to_load(end) == p.nt_to_recon_total
                if p.debug
                    fprintf('Finished reconstruting all %d volumes. \n', p.nt_to_recon);
                end
                break;
            end
            p.tpoints_to_load = p.tpoints_to_load(end) + (1:tpoints_step);
        else
            p.num_mux_cycle   = dat_num_mux_cycle;
            p.tpoints_to_load = p.skip_ref_pass + p.mux*p.num_mux_cycle + (1:tpoints_step);
        end
        
        if p.cal_flag
            p.cal_flag = 0;
            p.cal_dat_tpoints = [];
        end
                
    end  % -- End time point loop

    % -- Pseudo multiple replica, obtain noise map and retained SNR map
    if p.calc_snr_pmr
        pmr_rsnr = []; pmr_noise = []; pmr_snr = []; pmr_noise_ref = []; 
        if p.debug
            fprintf(' Pseudo-multiple replica simulation...\n');
            tpmr = tic;
        end
        if p.coil_compress
            p.num_coils = p.pfile_header.ncoils; % Change number of coils to before coil compression for pseudo-multiple replica simulation
        end
        [pmr_rsnr, pmr_noise, pmr_snr, pmr_noise_ref, nmsk] = simulate_pmr_snr(p, 1, p.pmr_raw_data); % pmr_rsnr Dim: [X(=p.nx_pres, if p.zpad_image=0; =p.zpad_size(1), if p.zpad_image=1), Y(=p.ny_pres, if p.zpad_image=0; =p.zpad_size(2), if p.zpad_image=1), Echo(=nec), Slice(=nz), Coil(=1), UndersamplingMask(=nmsk)].
        if p.debug
            tpmr = toc(tpmr);
            fprintf(' Pseudo-multiple replica simulation used %f s.\n', tpmr);
        end
        % - Rearrange the geometry of pmr results
        pmr_sl = cat(p.T_DIM, pmr_rsnr, pmr_noise);
        if p.debug
            pmr_sl = cat(p.T_DIM, pmr_sl, pmr_snr, pmr_noise_ref);
        end
        if bitand(p.dacq_ctrl, 4) % 3rd bit indicates odd-echo phase-flip (i.e., "pepolar")
            pmr_sl = flipdim(pmr_sl, p.PE_DIM);
        end
        if p.return_sos_im
            pmr_sl = reshape(pmr_sl, [nx_image, ny_image, p.mux, size(pmr_sl, p.T_DIM)]);
        else
            if strcmp(p.mux_recon_method, 'sense')
                pmr_sl = reshape(pmr_sl, [nx_image, ny_image, p.mux, size(pmr_sl, p.T_DIM)]);
            else
                if ~p.sense_coil_comb
                    pmr_sl = permute(reshape(pmr_sl, [nx_image, ny_image, p.mux, size(pmr_sl, p.C_DIM), size(pmr_sl, p.T_DIM)]), [1,2,3,5,4]); % Dim: [X, Y, Slice, Time, Coil]
                else
                    pmr_sl = reshape(pmr_sl, [nx_image, ny_image, p.mux, size(pmr_sl, p.T_DIM)]);
                end
            end
        end
        if p.transpose
            pmr_sl = permute(pmr_sl, [2,1,3,4,5]); % Dim: [1stImageDimension, 2ndImageDimension, Slice, Time, Coil(=nc_after_coil_compression if d contains complex images from GRAPPA recon; =1 otherwise)]
        end
        if p.rotation>0
            sz = [size(pmr_sl, 1), size(pmr_sl, 2), size(pmr_sl, 3), size(pmr_sl, 4), size(pmr_sl, 5)];
            if mod(p.rotation,2)
                d_size([1 2]) = d_size([2 1]);
                sz([1 2]) = sz([2 1]);
            end
            tmp = pmr_sl;
            pmr_sl = zeros(sz);
            for i=1:sz(3)
                for j=1:sz(4)
                    for k=1:sz(5)
                        pmr_sl(:,:,i,j,k) = rot90(tmp(:,:,i,j,k), p.rotation);
                    end
                end
            end
        end
        if p.swappf
            fprintf('Swapping freq/phase dims...\n');
            pmr_sl = permute(pmr_sl, [2,1,3,4,5]); % Dim: [1stImageDimension, 2ndImageDimension, Slice, Time, Coil(=nc_after_coil_compression if d contains complex images from GRAPPA recon; =1 otherwise)]
            d_size([1 2]) = d_size([2 1]);
        end

        % -- Scale the noise maps according to the receive gains
        pmr_sl(:, :, :, nmsk+1 : 2*nmsk, :) = pmr_sl(:, :, :, nmsk+1 : 2*nmsk, :) * scaling;
        if p.debug
            pmr_sl(:, :, :, end, :) = pmr_sl(:, :, :, end, :) * scaling;
        end
        
        % -- Save pmr results
        if save_vols && ~isempty(outfile)
            pmr_rsnr = pmr_sl(:, :, :, 1:nmsk, :);
            pmr_noise = pmr_sl(:, :, :, nmsk+1 : 2*nmsk, :);
            if p.debug
                pmr_snr = pmr_sl(:, :, :, 2*nmsk+1 : 3*nmsk, :);
                pmr_noise_ref = pmr_sl(:, :, :, end, :);
            end
            save(sprintf([outfile,'_pmr_s%03d'], sl-1), ...
                'pmr_rsnr', 'pmr_noise', 'pmr_snr', 'pmr_noise_ref', 'sl_loc', 'd_size', 'n_echo', 'ext_cal', 'pfile', '-v7');
        else
            if numel(p.slices_to_recon) > 1
                if isempty(pmr)
                    pmr = zeros([d_size(1:3), size(pmr_sl,4), size(pmr_sl,5)], 'like',pmr_sl);
                end
                pmr(:, :, sl_loc, :, :) = pmr_sl;
            else
                pmr = pmr_sl;
            end
        end
    end  % -- End If pseudo multiple replica
   
    if aborted_scan
        if p.debug
            fprintf('This is an aborted scan, reconstructed %d volumes. \n', p.nt_to_recon);
        end
    end

    % Put the slices into the right locations
    if ~save_vols && ~isempty(d_sl)
        if numel(p.slices_to_recon) > 1
            if isempty(d)
                d = zeros([d_size(1:3), size(d_sl,4), size(d_sl,5)], 'like',d_sl);
            end
            d(:, :, sl_loc, :, :) = d_sl;
        else
            d = d_sl;
        end
    end

end  % -- End slice loop

% Save all slices, all time points
if ~save_vols && (~isempty(d) && ~isempty(outfile))
    
    % Save calibration volumes in a seperate file
    d_tmp = d;
    [outpath, outname, ext] = fileparts(outfile);
    if p.internal_cal                                                   % save the first p.num_mux_cycle volumes for internal calibration
        d = d(:, :, :, 1:p.num_mux_cycle, :);
        outfile_cal = fullfile(outpath, [outname '_ical']);
    else
        switch p_cal.mux_excited
            case 1                                                      % save the first volume for single band external calibration
                d = d(:, :, :, 1, :);
                outfile_cal = fullfile(outpath, [outname '_sbcal']);
            case p.mux_excited                                          % save the first p_cal.num_mux_cycle volumes for external calibration
                d = d(:, :, :, 1:p_cal.num_mux_cycle, :);
                outfile_cal = fullfile(outpath, [outname '_ecal']);
        end
    end
    d_size(4) = size(d,4);
    save(outfile_cal, 'd', 'd_size', 'sl_loc', 'n_echo', 'cal_type', 'ext_cal', 'pfile', 'gitcommit', '-v7');    
            
    % Save accelerated scan volumes
    d = d_tmp(:, :, :, d_size(4)+1 : end, :);
    d_size(4) = size(d,4);
    save(outfile, 'd', 'd_size', 'sl_loc', 'n_echo', 'cal_type', 'ext_cal', 'pfile', 'gitcommit', '-v7');
    
    if p.calc_snr_pmr
        pmr_rsnr = pmr(:, :, :, 1:nmsk, :);
        pmr_noise = pmr(:, :, :, nmsk+1 : 2*nmsk, :);
        if p.debug
            pmr_snr = pmr(:, :, :, 2*nmsk+1 : 3*nmsk, :);
            pmr_noise_ref = pmr(:, :, :, end, :);
        end
        save([outfile, '_pmr'], 'pmr_rsnr', 'pmr_noise', 'pmr_snr', 'pmr_noise_ref', 'sl_loc', 'd_size', 'n_echo', 'ext_cal', 'pfile', '-v7');
    end
end

% -- Timing
if debug && ~isempty(outfile)
    time_used = toc(time_used);
    save([outfile,'_debug'], 'd', 'p', 'time_used', '-v7');
end

return
