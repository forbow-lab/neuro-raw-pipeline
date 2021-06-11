% Main function for unit testing of mux_epi_recon package.
clc;
clear;
close all;

% Reconstruction methods to test
meths.oneD_grappa = 1;
meths.slice_grappa = 1;
meths.split_slice_grappa = 1;
meths.oneD_grappa_sense_coil_comb = 1;
meths.slice_grappa_sense_coil_comb = 1;
meths.split_slice_grappa_sense_coil_comb = 1;
meths.hybrid_space_sense = 1;

% Reconstruction inputs
in.pfile = 'unit_test_data/ex10958_se5_P49664.7';
in.outfile = [];
in.ext_cal = [];
in.slices = [1,2];
in.nt_to_recon = 2;
in.n_vcoils = 10;
in.debug = 1;
in.apply_fermi = 0;
in.use_homodyne = 1;
in.notch_thresh = 0;

% Other params
main_fun = 'mux_epi_main'; % which mux_epi_recon main function to test, 'mux_epi_main_offline' or 'mux_epi_main_RT' or 'mux_epi_main'
prev_fname = ''; % .mat file name for previous test results
save_fname = ''; % .mat file name for saving test results
disp_res = true; % True: display reconed images

% Test
unit_test(meths, in, main_fun, prev_fname, save_fname, disp_res);
