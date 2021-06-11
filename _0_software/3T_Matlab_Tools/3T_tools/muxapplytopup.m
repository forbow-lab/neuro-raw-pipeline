function muxapplytopup(fileFwd, fieldmapfile)
%% muxapplytopup
%
% Perform "applytopup" susceptibility distortion correction on mux epi nii data
% from field maps generated using muxtopup
%
% muxapplytopup(fileFwd, fieldmapfile);
%
% Input:
%	fileFwd:       Name of input file containing PE forward images to correct (without .nii).
%	fieldmapfile:  Name of input file containing fieldmaps generated in muxtopup (without _fieldcoef.nii.gz)
% Output:
%	None explicitly. Creates one file called:
% 		{fileFwd}_fldcor.nii.gz
%
% REQUIRES FSL installed so 'topup' can be found in path...
%  eg: setenv('PATH',[getenv('PATH'),':/usr/local/fsl/bin'])
%

fsldir = getenv('FSLDIR');
if (exist(fsldir, 'dir') ~= 7) && (exist('/usr/local/fsl','dir') == 7),
	fsldir = '/usr/local/fsl';
	setenv('FSLDIR', fsldir);
	disp(['-- setting environment variable FSLDIR=',fsldir]);
end
setenv('FSLOUTPUTTYPE','NIFTI')

% Check for sufficient input arguments
tic;
if(nargin < 2), error('Require minimum of 2 input arguments!'), end

% Check existance of {fileFwd}.nii header to be corrected
if(~exist(strcat(fileFwd,'.nii'), 'file') == 2); 
  error('Nifti file %s with forward PE epi mux data not found.',strcat(fileFwd,'.nii'));
end

% Check existance of {fieldmapfile}_fieldcoef.nii.gz field map file from muxtopup
if(~exist(strcat(fieldmapfile,'_fieldcoef.nii.gz'), 'file') == 2); 
  error('Nifti file %s_fieldcoef.nii.gz with topup processed fieldmaps not found.',strcat(fieldmapfile,'.nii'));
end  

% Make system call to applytopup
cmd = sprintf('applytopup --imain=%s --inindex=1 --datain=%s_topup_param.txt --method=jac --topup=%s --out=%s_fldcor',fileFwd,fieldmapfile,fieldmapfile,fileFwd);
system(cmd);

% Report runtime and completion status
toc
disp(sprintf('MUXAPPLYTOPUP: "%s.nii" file processed to "%s_fldcor.nii.gz"',fileFwd,fileFwd));
