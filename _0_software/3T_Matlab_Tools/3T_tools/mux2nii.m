function d = mux2nii(matfile)
%% mux2nii
%
% Convert processed mux epi data .mat file into .nii file
%
% mux2nii(matfile);
%
% Input:
%	matfile : 	Name of .mat file from mux epi processing (with .mat)
%
% Output:
%	None explicitly. Creates a Nifti file called {pfile}.nii (ie.
%	P12345.7.nii)
%
% CHANGES:
%	20180515,Carl Helmick: modified Pfile naming convention to include date-time and Exam# preceeding Pnum.
%

[pfilePath pfileName pfileExt] = fileparts(matfile);

% Read in .mat file
  if(exist(matfile, 'file') == 2)
    load(matfile)
  else
    error('Matlab file %s with processed epi mux data not found.',matfile);
  end

% Read in pfile header using GE Orchestra package
  pfileHandle = GERecon('Pfile.Load', pfile);
  header = GERecon('Pfile.Header', pfileHandle);

% Write nii file
  d = d(:,:,:,3:end);   % Discard first 2 reference volumes (possibly only needed for internal ref scans)
  d = flipdim(d,2);     % Reorient matrix so correct nii orientation
  nii = make_nii(d);

% Amend template header to account for acquired scan
  nii.hdr.dime.pixdim(2)    = header.ImageData.pixsize_X;   % X resolution
  nii.hdr.dime.pixdim(3)    = header.ImageData.pixsize_Y;   % Y resolution
  nii.hdr.dime.pixdim(4)    = header.ImageData.slthick;     % Z resolution
  nii.hdr.dime.pixdim(5)    = header.ImageData.tr/1000000;  % TR resolution
  
  save_nii(nii,strcat(pfile,'.nii'));

% If a {pfile}_tensor.dat file exists, generate {pfile}.bvec {PFILE}.bval files
  tensorfile = strcat(pfile,'_tensor.dat');
  if(exist(tensorfile, 'file') == 2)
    bvec = importdata(tensorfile);
    ndiff = bvec(1);
    bvec = reshape(bvec(2:end),[3,ndiff]);
    bvec(1,:) = -bvec(1,:);
    nT2 = size(d,4) - ndiff;
    bval = 1000*[zeros(1,nT2), ones(1,ndiff)]; % B1000 assumed.. find pfile flag (user1?)
    bvec = [zeros(3,nT2), bvec];
    fid = fopen([pfile,'.bval'],'w');
    fprintf(fid,'%.0f ',bval);
    fprintf(fid,'\n');
    fclose(fid);
    fid = fopen([pfile,'.bvec'],'w');
    for i=1:3,
      fprintf(fid,'%.3f ',bvec(i,:));
      fprintf(fid,'\n');
    end
    fclose(fid);
  end

  return;
