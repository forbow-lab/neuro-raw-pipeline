function [ksp, p] = epi_pha_correct(ksp, p)
%
% function [ksp, pha_flt] = epi_pha_correct(ksp, pha_coe, p)
%
% Correct the 0-th and 1st order phases in the x-ky space for EPI data.
% Correspondingly, the constant phases of and the misalignment between the
% odd and even echoes in the k-space are corrected.
%
% Inputs
%   ksp     - Uncorrected k-space data. Dim: [Kx, Ky, Echo, Slice, Coil, Time]
%   pha_coe - Coefficients for x-ky phase correction. See get_pha_flt.m for details.
%             Dim: [2(0th order, 1st order), other dimensions(e.g. ny, nsl, nc].
%   p       - Parameter structure. See mux_epi_params.m for details. The following
%             fields are used in this function: FE_DIM, PE_DIM, EC_DIM, SL_DIM, C_DIM, T_DIM.
%
% Outputs
%   ksp     - Corrected k-space data (note that the correction is done in-place).
%   pha_flt - X-ky phase filter. xKy_corrected = xKy_uncorrected .* pha_flt.
%             Dim: [X, Ky, Echo(=1), Slice, Coil].
%
% (c) Kangrong Zhu,     Stanford University     July 2012

datsz = get_dat_sz(ksp, p);
if p.cal_flag
    p.pha_flt = get_pha_flt(p.pha_coe_default, datsz.x);                              % Dim: [X, Ky, Slice, Coil]
    p.pha_flt = reshape(p.pha_flt, [datsz.x, datsz.y, 1, datsz.sl, datsz.c]); % Dim: [X, Ky, Echo(=1), Slice, Coil]
end
for echo = 1 : datsz.ec
    for t = 1 : datsz.t
        ksp(:, :, echo, :, :, t) = fftc( ifftc(ksp(:, :, echo, :, :, t), p.FE_DIM) .* p.pha_flt(:, :, :, 1:size(ksp, 4), :), p.FE_DIM); % The same coefficients are used for every echo and every time point. When pha_flt is computed with single band reference data and ksp is data from accelerated time points, pha_flt and ksp may have different sizes in slice dimension, in which case only the first nslice of pha_flt is used.
    end
end

return
