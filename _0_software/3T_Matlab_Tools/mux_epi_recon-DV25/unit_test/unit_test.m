function unit_test(meths, in, main_fun, prev_fname, save_fname, disp_res)
%
% function unit_test(meths, in, main_fun, [prev_fname], [save_fname], [disp_res=false])
%
% Basic unit testing for mux_epi_recon package.
%
% Inputs
%   meths      - A structure holding boolean values indicating which reconstruction methods to test.
%                Must contain the following fields: oneD_grappa, slice_grappa, split_slice_grappa, 
%                oneD_grappa_sense_coil_comb, slice_grappa_sense_coil_comb, split_slice_grappa_sense_coil_comb, hybrid_space_sense.
%   in         - A structure holding all input parameters to the mux_epi_main function except the parameter for the reconstruction method.
%                See mux_epi_main_offline or mux_epi_main_RT for details.
%                Must contain the following fields: pfile, outfile, ext_cal, slices, nt_to_recon, n_vcoils, debug, apply_fermi, use_homodyne, notch_thresh.
%   main_fun   - Name of main function to test. 'mux_epi_main_offline' or 'mux_epi_main_RT' or 'mux_epi_main'.
%   prev_fname - Name of a .mat file containing reconstruction results from previous testing.
%   save_fname - Name of a .mat file to save reconstruction results from this testing into.
%   disp_res   - True: display reconstructed images and SNR-related variables.
%
% (c) Kangrong Zhu  Stanford University     Oct 2015

addpath('../');

if ~exist('disp_res', 'var') || isempty(disp_res)
    disp_res = false;
end
    
%% Method names & figure attributes
s = set_meth_strs();
if disp_res
    f = set_fig_attrib(meths, s);
else
    f = [];
end

%% Whether to compare with previous results
cmp_prev = (exist('prev_fname', 'var') && ~isempty(prev_fname) && (exist(prev_fname, 'file') || exist([prev_fname '.mat'], 'file')));
if cmp_prev
    load(prev_fname);
end

%% Whether to save this test's results
save_res = (exist('save_fname', 'var') && ~isempty(save_fname));
if save_res
    addpath('../utils/');
    time = get_current_time;
    test_start_time = [time{1} '-' time{2} '-' time{3} ' ' time{4} ':' time{5} ':' time{6}];
    save(save_fname, 'test_start_time');
end

%% Test method by method
for meth_idx = 1 : s.tot_num_meths
    meth = s.varnames{meth_idx};
    
    % Rename previous results
    cmd = sprintf([...
        'if cmp_prev && exist(''d_' meth ''', ''var'')\n',...
        '    d_' meth '_prev = d_' meth ';\n',...
        'else\n',...
        '    d_' meth '_prev = [];\n',...
        'end\n']);
    eval(cmd);

    cmd = sprintf([...
        'if cmp_prev && exist(''snr_d_' meth ''', ''var'')\n',...
        '    snr_d_' meth '_prev = snr_d_' meth ';\n',...
        'else\n',...
        '    snr_d_' meth '_prev = [];\n',...
        'end\n']);
    eval(cmd);

    % Test
    cmd = sprintf([...
        'if meths.' meth '\n',...
        '    in.recon_method = ''' s.varname2input(meth) ''';\n',...
        '    [d_' meth ', snr_d_' meth ', f] = test_one_method(main_fun, in, d_' meth '_prev, snr_d_' meth '_prev, f);\n\n',...
        '    if save_res\n',...
        '        save(save_fname, ''d_' meth ''', ''snr_d_' meth ''', ''-append'');\n',...
        '    end\n',...
        'end\n']);
    eval(cmd);
end

return

function [d, snr_d, fig] = test_one_method(main_fun, in, d_prev, snr_d_prev, fig)
%
% Test one reconstruction method.
%
% Inputs
%   main_fun   - Main function name to test. See inputs to function unit_test for details.
%   in         - A structure holding all inputs to the main function. See inputs to function unit_test for details.
%   d_prev     - Previously reconstructed images.
%   snr_d_prev - Previously calculate SNR attributes.
%   fig        - A structure holding figure attributes. Output of function set_fig_attrib (display) or empty matrix (no display).

%% Reconstruct
switch main_fun
    case 'mux_epi_main'
        snr_out_str = ', snr_d';
    case 'mux_epi_main_offline'
        snr_out_str = ', snr_d';
    case 'mux_epi_main_RT'
        snr_out_str = '';
end
cmd = ['[d' snr_out_str '] = ' main_fun '(in.pfile, in.outfile, in.ext_cal, in.slices, in.nt_to_recon, in.n_vcoils, in.debug, in.recon_method, in.apply_fermi, in.use_homodyne, in.notch_thresh);'];
eval(cmd);

%% Display
if ~isempty(fig)
    d = double(d);
    if ~exist('snr_d', 'var')
        snr_d = [];
    end
    fig.plots = fig.plots + 1;

    % String for reconstruction method
    D_PREV_INPUT_ORDER = 3;
    recon_method_str = inputname(D_PREV_INPUT_ORDER);
    recon_method_str = strsplit(recon_method_str, '_');
    recon_method_str(1) = []; % Deletes first cell, which is 'd'
    recon_method_str(end) = []; % Deletes last cell, which is 'prev'
    recon_method_str = strjoin(recon_method_str, ' '); % String to be used for variable names
    
    % Variables to display
    vars{1} = 'd';
    if ~isfield(fig.dispr, 'im') || isempty(fig.dispr.im)
        fig.dispr.im = [0, fig.dispr.im_max_ratio * max(abs(d(:)))];
    end
    dispr{1} = fig.dispr.im;
    cmap{1} = fig.cmap.im;
    
    if ~isempty(d_prev) && (length(size(d)) == length(size(d_prev))) && all(size(d) == size(d_prev))
        d_difference = d - d_prev;
        vars{end+1} = 'd_difference';
        dispr{end+1} = fig.dispr.im .* fig.diff_wd_ratio;
        cmap{end+1} = fig.cmap.im;
        cmd = ['fprintf([''' recon_method_str ' reconed: max abs difference between current and previous test is %.3f\n''], max(abs(d_difference(:))));'];
        eval(cmd);
    end
    
    if exist('snr_d', 'var') && ~isempty(snr_d)
        snr_d_fields = {'sense_rsnr', 'pmr_rsnr', 'pmr_snr', 'amr_snr'};
        dispr_ranges = {fig.dispr.rsnr, fig.dispr.rsnr, fig.dispr.snr, fig.dispr.snr};
        cmap_types = {fig.cmap.rsnr, fig.cmap.rsnr, fig.cmap.snr, fig.cmap.snr};
        
        for fidx = 1 : length(snr_d_fields)
            field_name = snr_d_fields{fidx};
            cmd = sprintf([...
                'if isfield(snr_d, ''' field_name ''') && ~isempty(snr_d.' field_name ')\n',...
                '    vars{end+1} = ''snr_d.' field_name ''';\n',...
                '    dispr{end+1} = dispr_ranges{fidx};\n',...
                '    cmap{end+1} = cmap_types{fidx};\n\n',...
                '    if ~isempty(snr_d_prev) && isfield(snr_d_prev, ''' field_name ''') && ~isempty(snr_d_prev.' field_name ') && (length(size(snr_d.' field_name ')) == length(size(snr_d_prev.' field_name '))) && all(size(snr_d.' field_name ') == size(snr_d_prev.' field_name '))\n',...
                '        vars{end+1} = ''snr_d.' field_name ' - snr_d_prev.' field_name ''';\n',...
                '        dispr{end+1} = dispr_ranges{fidx} .* fig.diff_wd_ratio;\n',...
                '        cmap{end+1} = cmap_types{fidx};\n',...
                '    end\n',...
                'end\n']);
            eval(cmd);                
        end
    end
    
    % Display variables
    for vidx = 1 : length(vars)
        var_name = vars{vidx};
        figure(fig.num(var_name));
        subplot(fig.nrows, fig.ncols, fig.plots);
        cmd = ['imshowALL(abs(' var_name ') + eps, dispr{vidx} , [1, size(' var_name ', 3)*size(' var_name ', 4)]);'];
        eval(cmd);
        colormap(cmap{vidx});
        tit = title([recon_method_str, ': ' var_name]);
        set(tit, 'interpreter', 'none');
    end
end

return

function s = set_meth_strs()
% Set up a structure holding strings for reconstruction method names.

% Strings for variable names
s.varnames = {'oneD_grappa',...
    'slice_grappa',...
    'split_slice_grappa',...
    'oneD_grappa_sense_coil_comb',...
    'slice_grappa_sense_coil_comb',...
    'split_slice_grappa_sense_coil_comb',...
    'hybrid_space_sense'};

% Strings for inputs into reconstruction package
input_meth_names = {'1Dgrappa',...
    'slice-grappa',...
    'split-slice-grappa',...
    '1Dgrappa_sense1',...
    'slice-grappa_sense1',...
    'split-slice-grappa_sense1',...
    'sense'};

s.varname2input = containers.Map(s.varnames, input_meth_names); % containers.Map(keySet, valueSet);

% Total number of possible methods
s.tot_num_meths = s.varname2input.Count;

return

function f = set_fig_attrib(meths, s)
% Set up a structure holding attributes of the figure for displaying results.

keySet = {'d',...
    'd_difference',...
    'snr_d.sense_rsnr',...
    'snr_d.pmr_rsnr',...
    'snr_d.pmr_snr',...
    'snr_d.amr_snr',...
    'snr_d.sense_rsnr - snr_d_prev.sense_rsnr',...
    'snr_d.pmr_rsnr - snr_d_prev.pmr_rsnr',...
    'snr_d.pmr_snr - snr_d_prev.pmr_snr',...
    'snr_d.amr_snr - snr_d_prev.amr_snr'};
valueSet = 500 + (1 : length(keySet));
f.num = containers.Map(keySet, valueSet);

f.nrows = 0;
for meth_idx = 1 : s.tot_num_meths
    meth = s.varnames{meth_idx};
    cmd = sprintf([...
        'if meths.' meth '\n', ...
        '    f.nrows = f.nrows + 1;\n', ...
        'end']);
    eval(cmd);
end

f.ncols = 1;
f.plots = 0;

f.dispr.im_max_ratio = 0.6;
f.dispr.rsnr = [0, 1.5];
f.dispr.snr = [0, 100];
f.cmap.im = gray(256);
f.cmap.rsnr = jet(256);
f.cmap.snr = jet(256);

f.diff_wd_ratio = 0.01;

return