function muxrecon(pfile, outfile, ext_cal, slice, n_vcoils, recon_method, vol_flag, GPU_flag)

if(~exist('vol_flag', 'var'));    vol_flag = '1';     end
if(~exist('GPU_flag', 'var'));    GPU_flag = '1';     end

fprintf('mux_epi_main_RT("%s", "%s", "%s", %d, [], %d, 0, "%s", 1, 1, 0, %d, %d)\n', pfile, outfile, ext_cal, str2num(slice), ...
    str2num(n_vcoils), recon_method, str2num(vol_flag), str2num(GPU_flag));
mux_epi_main_RT(pfile, outfile, ext_cal, str2num(slice), [], str2num(n_vcoils), 0, recon_method, 1, 1, 0, str2num(vol_flag), str2num(GPU_flag));    


