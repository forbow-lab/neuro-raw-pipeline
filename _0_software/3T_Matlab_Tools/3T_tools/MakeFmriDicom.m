function MakeFmriDicom(t_thresh, series_num, series_name, dir_anat_dicom, file_func_nii, dir_out_dicom)
%MakeFmriDicom Create dicom with merged anatomic and function for surgical planning
%   MakeFmriDicom(dir_anat_dicom, file_func_nii, dir_out_dicom);
%
%   Input:
%      t_thresh:       Threshold value for t-stat functional map to use
%      series_num:     Series number for fused dicom series
%      series_name:    Series name for fused dicom series
%      dir_anat_dicom: Directory containing T1w anatomical dicom files (Default queries folder)
%      file_func_nii:  File containing functional map .nii file (Default queries folder)
%      dir_out_dicom:  Directory containing T1w anat with functional overlay (Default queries name)
%   Output:
%      None Explicitly.
%      Dicom folder containing fused T1w anat (thresh 0-2000) and functional overlay (thresh 20000-30000) created in local directory

DEBUG = 1;  % Debug orientation plot

% Arguments check
if (nargin < 5)
   disp('No t-threshold, series_num, dicom and nifti folder names provided. Prompting for input.')
   
   default_out = num2str('2.3');
   prompt={'Threshold value for t-stat map? :'};
   dlg_title='Input';
   num_lines=1;
   def={default_out};
   t_thresh = inputdlg(prompt,dlg_title,num_lines,def);
   t_thresh = str2num(t_thresh{1});
   fprintf('T-Threshold selected: %4.2f\n',t_thresh)
   
   default_out = num2str('100');
   prompt={'Series number for fused series? :'};
   dlg_title='Series Number';
   num_lines=1;
   def={default_out};
   series_num = inputdlg(prompt,dlg_title,num_lines,def);
   series_num = str2num(series_num{1});
   fprintf('Series Number: %4.0d\n',series_num)
   
   default_out = num2str('Function Map');
   prompt={'Series name for fused series? :'};
   dlg_title='Series Name';
   num_lines=1;
   def={default_out};
   series_name = inputdlg(prompt,dlg_title,num_lines,def);
   series_name = series_name{1};
   fprintf('Series Name: %s\n',series_name)
   
   dir_anat_dicom = uigetdir(pwd,'Folder containing anatomic dicom images');
   fprintf('Input Anatomic Dicom Dir: %s\n',dir_anat_dicom);

   file_func_nii = uigetfile('*.nii','Nifti file containing functional threshold map');
   fprintf('Input Functional Nifti Map Filename: %s\n',file_func_nii);

   default_out = strcat(dir_anat_dicom,'_',strrep(series_name,' ','_'));
   prompt={'Output dicom directory name :'};
   dlg_title='Input';
   num_lines=1;
   def={default_out};
   dir_out_dicom = inputdlg(prompt,dlg_title,num_lines,def);
   dir_out_dicom = dir_out_dicom{1};
   dir_out_dicom = strrep(dir_out_dicom,' ','_');
   fprintf('Output Fused Dicom Folder Name: %s\n',dir_out_dicom)
else
   disp('Using input values provided below.');
   fprintf('T-Threshold selected: %4.2f\n',t_thresh)
   fprintf('Series Number: %4.0d\n',series_num)
   fprintf('Series Name: %s\n',series_name)
   fprintf('Input Anatomic Dicom: %s\n',dir_anat_dicom);
   fprintf('Input Functional Nifti Map: %s\n',file_func_nii);
   fprintf('Output Fused Dicom Folder Name: %s\n',dir_out_dicom)
end

% Load functional nifti volume
data_func = load_nii(file_func_nii);
data_func = data_func.img;

% Create anatomic nifti volume
tmpimdir = [pwd,'/tmpimdir/'];
if(exist(tmpimdir,'dir'))
    rmdir(tmpimdir,'s');
end
mkdir(tmpimdir);
CMD = sprintf('/biotic/opt/mricron/dcm2nii -o "%s" -b /biotic/home/chrisb/matlab/ACR_Analysis/dcm2nii.ini "%s"', ...
    tmpimdir,dir_anat_dicom);
[status, result] = system(CMD);
if(status)
    disp(result)
    error('Problem running dcm2nii dicom to nifti conversion');
end
files = dir(tmpimdir);
files = files(~[files.isdir]);  % Select only real files (not hidden or sub directories)
filenii = [tmpimdir files(1).name];
data_anat = load_nii(filenii);
data_anat = data_anat.img;

% Read dicom slice array
dicomdict('set','/lauterbur/opt/ESE/ESE_DV25.0_R02/DICOM/examples/matlab/gems-dicom-dict.txt');
fn = dir(fullfile(dir_anat_dicom, 'MR*'));

for f = 1:length(fn)
  file = fullfile(dir_anat_dicom,fn(f).name);
  data_anat_dcm(:,:,f) = dicomread(file);
  info = dicominfo(file,'UseDictionaryVr',true);
end

% Scale anatomic and functional maps to fit into 0-2000, 20000-30000 range
data_anat = double(data_anat);
data_func = double(data_func);
data_anat_sc = data_anat*2000.0/max(data_anat(:));
data_func_sc = (data_func-t_thresh)*10000/(max(data_func(:))-t_thresh) + 20000;
ind = data_func > t_thresh;
data_fused = data_anat_sc;
data_fused(ind) = data_func_sc(ind);
data_fused = int16(data_fused);

% Adjust fused images for dicom orientation
data_fused = flipud(fliplr(permute(data_fused,[2 1 3])));

% Plot optional figure to ensure orientation correct
if(DEBUG)
    figure;
    colormap('gray')
    maxsig = 0.8*max(data_anat_dcm(:));
    subplot(2,2,1);
    imagesc(data_anat_dcm(:,:,floor(0.25*length(fn))),[0 maxsig])
    subplot(2,2,2);
    imagesc(data_fused(:,:,floor(0.25*length(fn))),[0 0.8*2000])
    subplot(2,2,3);
    imagesc(data_anat_dcm(:,:,floor(0.75*length(fn))),[0 maxsig])
    subplot(2,2,4);
    imagesc(data_fused(:,:,floor(0.75*length(fn))),[0 0.8*2000])
end

% Read dicom array and generate fused output
NewSeriesUID = dicomuid;   % Is this OK?
if(exist(dir_out_dicom))
    rmdir(dir_out_dicom,'s');
end
mkdir(dir_out_dicom);
warning('off','all');
for f = 1:length(fn)
  file_in = fullfile(dir_anat_dicom,fn(f).name);
  file_out = fullfile(dir_out_dicom,fn(f).name);
  info = dicominfo(file_in,'UseDictionaryVr',true);
  info.SeriesDescription = series_name;
  info.SeriesNumber = series_num;
  info.SeriesInstanceUID = NewSeriesUID;
  dicomwrite(data_fused(:,:,f), file_out, info, 'WritePrivate', true);
end
warning('on','all');
dicomdict('factory');
rmdir(tmpimdir,'s');

% Ask if dicoms should be pushed to AW Server?
default_out = num2str('y');
prompt={'Push dicom to AW Server (y/n)? :'};
dlg_title='(y/n)?';
num_lines=1;
def={default_out};
pushdcm = inputdlg(prompt,dlg_title,num_lines,def);
pushdcm = pushdcm{1};
if (strcmp(pushdcm(1),'y'))
    pushDICOMs(dir_out_dicom);
end
