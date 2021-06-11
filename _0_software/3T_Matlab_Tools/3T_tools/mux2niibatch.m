function mux2niibatch(pfiledir)
%% mux2niibatch
%
% Convert all mux epi P-file data found in "path" to .nii files
%
% mux2niibatch([path]);
%
% Input:
%	path : 	Name of directory for processing all mux epi p-files to nifti
%           (Default is current working directory)
%
% Output:
%	None explicitly. Creates Nifti files called {pfile}.nii for all p-files in path
%  (ie.	20180515-1836_E01234_P12345.7.nii)
%
% REQUIRES FSL installed and added to PATH variable
%  eg: setenv('PATH',[getenv('PATH'),':/usr/local/fsl/bin'])
%
% CHANGES:
%	20180515,Carl Helmick: modified Pfile naming convention to include date-time and Exam# preceeding Pnum.
%

fsldir = getenv('FSLDIR');
if (exist(fsldir, 'dir') ~= 7) && (exist('/usr/local/fsl','dir') == 7),
	fsldir = '/usr/local/fsl';
	setenv('FSLDIR', fsldir);
	disp(['-- setting environment variable FSLDIR=',fsldir]);
end
if (exist(fsldir,'dir') == 7),
    %disp(['fsldir=',fsldir]);
    fsldirmpath = sprintf('%s/etc/matlab',fsldir);
    %disp(['fsldirmpath=',fsldirmpath]);
    path(path,fsldirmpath);
    setenv('FSLOUTPUTTYPE','NIFTI');
    setenv('PATH',[getenv('PATH'),':',fsldir,'/bin']);
    clear fsldir fsldirmpath;
else
    disp(['***ERROR: cannot locate environment FSLDIR or /usr/local/fsl...exiting']);
    exit(1);
end


% Locate all p-files within specified "pfiledir"
if(~exist('pfiledir', 'var')),
  pfiledir = pwd;
end
pfiles = dir(strcat(pfiledir,'/20*_E*_P*.7'));	%pattern: "20180509-1336_E09625_P44544.7"
numfiles = length(pfiles);
if(numfiles == 0), 
    disp(sprintf('***ERROR: found no P-files in specified directory = ''%s''',pfiledir)); 
    return;
end

% Loop all p-files found and process into .mat files
  for ind = 1:numfiles,
      pfile = pfiles(ind).name;
      if exist([pfile,'.nii'], 'file') == 2 || exist([pfile,'.nii.gz'], 'file') == 2,
      	disp(['-- found existing NIFTI file for pfile=',pfile,'...skipping.']);
      	continue
      end
      disp(['-- converting pfile: ',pfile]);
      try,
     	 mux_epi_main_offline(pfile,strcat(pfile,'.mat'));
      catch,
      	fprintf('WARNING: caught error attempting to run mux_epi_main_offline(pfile,pfile.mat), skipping pfile=%s\n\n', pfile);
      	continue
      end
      mux2nii(strcat(pfile,'.mat'));
  end
  
return
