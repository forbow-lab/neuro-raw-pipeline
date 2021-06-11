function muxtopup(fileFwd, fileRev, outfile, ncal)
%% muxtopup
%
% Perform "topup" susceptibility distortion correction to create field maps on mux epi nii data
%
% muxtopup(fileFwd, fileRev, outfile, ncal);
%
% Input:
%	fileFwd:   Name of input file containing PE forward images (without .nii).
%	fileRev:   Name of input file containing PE reversed images (without .nii)
%	outfile:   Base name of field map coefficient and motion files. Default "fieldmap"
%	ncal:      Use 'ncal' first time points (if present) from fileFwd and fileRev
%		       (Default uses ncal=4, or all available if fewer)
% Output:
%	None explicitly. Creates three files called:
% 		{outfile}_fieldcoef.nii.gz
%		{outfile}_movepar.txt
%       {outfile}_topup_param.txt
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

DEFAULT_NCAL=4;

% Check for sufficient input arguments
  tic
  if(nargin < 2), error('Require minimum of 2 input arguments!'), end
  
% Read in {fileFwd}.nii header
  if(exist(strcat(fileFwd,'.nii'), 'file') == 2)
    fhdr = load_nii_hdr(strcat(fileFwd,'.nii'));
  else
    error('Nifti file %s with forward PE epi mux data not found.',strcat(fileFwd,'.nii'));
  end
  
% Read in {fileRev}.nii header
  if(exist(strcat(fileRev,'.nii'), 'file') == 2)
    rhdr = load_nii_hdr(strcat(fileRev,'.nii'));
  else
    error('Nifti file %s with reverse PE epi mux data not found.',strcat(fileRev,'.nii'));
  end
  
% Check if optional input ncal and outfile exist, otherwise set to default
  if(~exist('ncal','var'));	ncal = DEFAULT_NCAL;		end
  if(~exist('outfile','var'));	outfile = 'fieldmap';	end

% Set ncal to minimum number of volumes from fileFwd, fileRev or to ncal
  if(fhdr.dime.dim(5) < ncal), ncal = fhdr.dime.dim(5);	end
  if(rhdr.dime.dim(5) < ncal), ncal = rhdr.dime.dim(5);	end

% Check that even number of slices for topup, or truncate last slice
  if(mod(fhdr.dime.dim(4),2)), 
      disp('Warning. Odd number of slices, last trimmed for topup processing.');
      nsf = fhdr.dime.dim(4)-1; end
  if(mod(rhdr.dime.dim(4),2)), nsr = fhdr.dime.dim(4)-1; end
  nx = fhdr.dime.dim(2); ny = fhdr.dime.dim(3);

% Truncate fimg and rimg to ncal dimension
  system(sprintf('fslroi %s fimg 0 %d 0 %d 0 %d 0 %d',fileFwd,nx,ny,nsf,ncal));
  system(sprintf('fslroi %s rimg 0 %d 0 %d 0 %d 0 %d',fileRev,nx,ny,nsr,ncal));
  system('fslmerge -t bimg fimg rimg');

% Create topup_param.txt file to control field map (topup) calculation
  [fid,status] = fopen(sprintf('%s_topup_param.txt',outfile),'w');
  if(fid<0), error('Failed to open %s_topup_param.txt due to message: %s',outfile,status); end
  for ind = 1:ncal,
    fprintf(fid,'0 1 0 1\n');
  end
  for ind = 1:ncal,
    fprintf(fid,'0 -1 0 1\n');
  end
  fclose(fid);

% Make system call to topup
  cmd = sprintf('topup --imain=bimg --datain=%s_topup_param.txt --config=b02b0.cnf --out=%s',outfile,outfile);
  system(cmd);
  
% Clean up intermediate processing files and report usage time and completion
  system('/bin/rm -rf fimg.nii.gz rimg.nii.gz bimg.nii.gz bing.topup_log');
  toc
  disp(sprintf('MUXTOPUP: file "%s" processed',outfile));

