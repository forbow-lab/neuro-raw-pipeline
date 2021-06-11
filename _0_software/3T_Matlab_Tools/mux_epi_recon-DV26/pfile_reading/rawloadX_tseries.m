function [data, params, header] = rawloadX_tseries(fname,frames,echoes,slices,coils,passes,params)
%
% Modified from rawloadX.m, to load k-space data of a time series conntained
% in a p-file.
%
% function [data,header] = rawloadX_tseries(fname,frames,echoes,slices,coils,passes[,params]);
%
%	Function reads selected (or all) data from a P-file,
%	using given lists of frames, echoes, slices and/or coils.
%	Reads 12.x, 14.x and 20.x EPIC files.
%	See 'help geX' for more information.
%
%	INPUT:
%		fname   = path and file name.
%		frames  = frame numbers to read (+hnover) (Frame 0=baseline).
%		echoes  = echo numbers to read (1...).
%		slices  = slice numbers to read (1...).
%		coils   = coil numbers to read (1...).
%               passes  = pass numbers to read (1...). For the fMRI scan,
%                         this must be continuous from 1 to the last pass
%                         number to read. For an external calibration scan,
%                         this should correspond to the last group of slice
%                         phase cycling.
%               params  = The structure with all parameters of the scan and
%                         for the reconstruction. Please refer to mux_epi_params.m
%                         for the definition of each field in this structure.
%                         The 'pfile_header' field in 'params' is used in
%                         this function. If 'params.pfile_header' is empty, the
%                         header information will be read using function 'rawheadX.m'.
%
%	OUTPUT:
%		data = data array (up to 6 dimensions, including all
%			specified frames, echoes, slices, coils, passes.
%               params_out = The output parameter structure, having the same
%                       fields as the input 'params'. The following fields might
%                       have been changed inside this function if the scan was aborted
%                       midway and the p-file was only partially collected: num_mux_cycle,
%                       num_passes, nt_to_recon, tpoints_to_load.

%		header = header structure - see rawheadX.m.
%
%	EXAMPLES:
%	  d = rawloadX('P00000.7');	% Read full p-file, no baselines.
%	  d = rawloadX('P00000.7',[],[],[],1);	% Read only 1st coil.
%	  d = rawloadX('P00000.7',[0:1024]);	% Read all data with baselines
%						% (like rawload, rawload12X etc)
%	  d = rawloadX('P00000.7',[],[],[1,2,4,7]; % Read slices 1,2,4,7.
%
%	NOTES:
%	1) If ranges are beyond those of the header, then they are
%		limited by ranges in the header.
%	2) If not all expected frames are read, then the data are returned
%		as a 2D array, with all frames read.
%	3) See rawheadX to just read header information.
%
%	SEE ALSO:  rawheadX, writepfileX, reconX
%
%	B.Hargreaves -- April 2008.
%       Modified by Kangrong Zhu,   Stanford University       July 2012



% -- Set defaults, if not provided.
if (nargin < 2) frames=[]; end;
if (nargin < 3) echoes=[]; end;
if (nargin < 4) slices=[]; end;
if (nargin < 5) coils=[];  end;
if (nargin < 6) passes=[]; end;

% -- Header information
if isempty(params.pfile_header)
    header = rawheadX(fname);
else
    header = params.pfile_header;
end
if isfield(params, 'precision'), precision = params.precision; else precision = 'single'; end
if isfield(params, 'slice_shuffling'), slice_shuffling = params.slice_shuffling; else slice_shuffling = 0; end
if isfield(params, 'mux'), mux = params.mux; else mux = 1; end
if isfield(params, 'num_mux_cycle'), num_mux_cycle = params.num_mux_cycle; else num_mux_cycle = 0; end

% -- Open file.
fip = fopen(fname,'r','l');
if fip == -1
  error('File %s not found\n',fn);
end

% -- Check parameters passed are reasonable.
if (isempty(frames)) frames = 1 : header.nframes+header.hnover;  end;
if (isempty(echoes)) echoes = 1 : header.nechoes;                end;
if (isempty(slices)) slices = 1 : header.nslices/header.npasses; end;
if (isempty(coils))  coils  = 1 : header.ncoils;                 end;
if (isempty(passes)) passes = 1 : header.npasses;                end;

frames = checkrange(frames,0,header.nframes+header.hnover,'Frames');
echoes = checkrange(echoes,1,header.nechoes,'Echoes');
slices = checkrange(slices,1,header.nslices/header.npasses,'Slices');
coils  = checkrange(coils,1,header.ncoils,'Coils');
passes = checkrange(passes,1,header.npasses,'Passes');

%if (passes(1)==1) && (passes(end) ~= length(passes)) % Assuming num_mux_cycle is always larger than 1, so the p-file must be from an fMRI scan, not from an external calibration scan, when passes(1)==1
%    error('The pass numbers to read must be continuous from 1 to the last pass number to read.');
%end

% --- Allocate array for data, based on passed arguments or defaults.
data = zeros(header.frsize,length(frames),length(echoes),length(slices),length(coils),length(passes), precision);

% --- Change all to start at 0, EXCEPT frames, as the default is to ignore the baseline.
echoes = echoes-1;
slices = slices-1;
coils  = coils-1;
passes = passes-1;

ptsize = header.ptsize;			% Sample size (bytes)
framesize = 2*ptsize*header.frsize;	% Frame size in bytes.
readsize = framesize/ptsize;
echosize = framesize*(1+header.nframes+header.hnover); % Size of one echo (bytes)
slicesize = echosize*header.nechoes;	% Size of one slice (bytes)
coilsize = slicesize*header.nslices/header.npasses; % Size of one echo (bytes)
passsize = coilsize*header.ncoils;      % Size of one pass (bytes)
dattypes = {'int16','int32'};
datatype = dattypes{ptsize/2};

rawhdrsize = header.rawhdrsize;
rawdatasize = header.rawsize;

aborted_scan = false;                   % true: the scan was aborted midway and the k-space data was collected only partially.

% skip past the header (to the DISDAQ frame)
status = fseek(fip,rawhdrsize,-1);

% --- For each entity (pass, coil, slice, echo, frame) we have variables
%	count and index.  Index indicates the next index into the
%	list for the entity (ie next coil is coils(coilindex).
%	count indicates the file position for that entity.

totalframes = 0;	% Desired frames
framespresent=0;	% Frames successfully read.

passcount = 0;
passindex = 1;
while(passindex <= length(passes))
    passskip = passes(passindex)-passcount;
    status = fseek(fip,passskip*passsize,0); % Skip to desired pass
    if status~=0 || feof(fip)
        aborted_scan = true;
        break;
    end
    passcount = passcount+passskip;

    coilcount = 0;
    coilindex = 1;		% Index into coils.
    while (coilindex <= length(coils))
        [fip,coilcount,aborted_scan] = skipsection(fip,coilcount,coils(coilindex),coilsize);
        if aborted_scan, break; end

        sliceindex = 1;		% Index into slices
        slicecount = 0;
        while (sliceindex <= length(slices))
            slices_act = slices;
            if (slice_shuffling == 1) 
                if params.isdifscan 
                    if (passes(passindex) >= params.mux*params.num_mux_cycle) 
                        shift = passes(passindex) - params.mux*params.num_mux_cycle;
                    else
                        shift = 0;
                    end
                else
                    shift = passes(passindex);
                end
                slices_act = mod(slices - shift, params.num_slices);
            end

            [fip,slicecount,aborted_scan] = skipsection(fip,slicecount,slices_act(sliceindex),slicesize);
            if aborted_scan, break; end

            echoindex = 1;		% Index into echoes
            echocount = 0;
            while (echoindex <= length(echoes)) % -- Skip echoes before desired echo.
                [fip,echocount,aborted_scan] = skipsection(fip,echocount, echoes(echoindex),echosize);
                if aborted_scan, break; end

                frameindex = 1;		% Index into frames
                framecount = 0;
                while (frameindex <= length(frames))
                    % -- Skip frames before desired frame.
                    [fip,framecount,aborted_scan] = skipsection(fip,framecount,frames(frameindex),framesize);
                    if aborted_scan, break; end

                    % DEBUG:tt = sprintf('Reading Frame %d, Echo %d, Slice %d, Coil %d, Pass %d', frames(frameindex),echoes(echoindex),slices(sliceindex),coils(coilindex),passes(passindex)); disp(tt);

                    % Read frame
                    [dr,nread] = fread(fip,readsize,datatype);
                    totalframes = totalframes+1;	% Expected frames.

                    if (nread == readsize) && (status == 0)
                        dr = complex(dr(1:2:end),dr(2:2:end));
                        framespresent = framespresent+1;	% Keep track.
                    else
                        aborted_scan = true;
                        acquired_num_passes = passindex - 1;    % Assuming the data untill the previous pass has been collected correctly.
                        break;                                  % Once we know the p-file was only collected partially, we don't need to read in any more data.
                    end;

                    % -- Place in data array.
                    data(:,frameindex,echoindex,sliceindex,coilindex,passindex) = dr;  % -- Convert to complex.

                    framecount=framecount+1; frameindex=frameindex+1;
                end;  % --Frame loop.

                if aborted_scan
                    break;   % Once we know the p-file was only collected partially, we don't need to read in any more data.
                end

                % -- Skip remaining frames.
                [fip,framecount,aborted_scan] = skipsection(fip,framecount, header.nframes+header.hnover+1,framesize);
                if aborted_scan, break; end

                echocount=echocount+1; echoindex=echoindex+1;
            end;  % --Echo loop.

            if aborted_scan
                break;
            end

            % -- Skip remaining echoes
            [fip,echocount,aborted_scan] = skipsection(fip,echocount, header.nechoes, echosize);
            if aborted_scan, break; end

            slicecount=slicecount+1; sliceindex=sliceindex+1;

        end;  % --Slice loop.

        if aborted_scan
            break;
        end

        % -- Skip remaining slices
        [fip,slicecount,aborted_scan] = skipsection(fip,slicecount,header.nslices/header.npasses,slicesize);
        if aborted_scan, break; end

        coilcount=coilcount+1; coilindex=coilindex+1;
    end;  % --Coil loop.

    if aborted_scan
        break;
    end

    % -- Skip remaining coils
    [fip,coilcount,aborted_scan] = skipsection(fip,coilcount,header.ncoils,coilsize);
    if aborted_scan, break; end

    passcount=passcount+1; passindex=passindex+1;
end % --Pass loop
% -- No need to skip remaining passes, as no more data to read!

fclose(fip);	% -- Close file.

if ~aborted_scan
    acquired_num_passes = passindex - 1;
    params.tpoints_to_load = passes(1) + (1 : acquired_num_passes);
else
    acquired_num_passes = 0;
end

% -- Deal with partially acquired p-files.
if isfield(params, 'cal_flag') && params.cal_flag && acquired_num_passes < 1          % p-file contains no data.
    warning('The p-file doesn''t contain any data.');
end

if isfield(params, 'cal_flag') && params.cal_flag && acquired_num_passes < mux*num_mux_cycle   % Data in p-file are all mux phase cycling frames.
    warning(['The number of temporal frames contained in the p-file is %d, ', ...
        'mux is %d. Can''t conduct any reconstruction in this case.'], acquired_num_passes, mux);
end

% params.num_passes = acquired_num_passes;
params.nt_to_recon = acquired_num_passes + passes(1) - mux*num_mux_cycle;
if acquired_num_passes < length(passes)
    data = data(:, :, :, :, :, 1:acquired_num_passes);
end

if params.debug
    fprintf('  Only the acquired %d frames will be loaded and reconstructed.\n', params.num_passes);
end

return;		% -- End of main function.



function arrout = checkrange(arr,amin,amax,atype)
%%%% -- Check values in array are within range, and remove ones that
	%%%are not!
f = find(arr >= amin);
arrout = arr(f);
f = find(arrout <= amax);
arrout = arrout(f);

if (length(arrout) < length(arr))
  arr = arr(:);
  f1 = find(arr < amin);
  f2 = find(arr > amax);
  tt = sprintf('%d %s out of range: ',length(f1)+length(f2),atype); disp(tt);
  disp([arr(f1); arr(f2)].');
end;
return;
%%%% --------


function [fileptr,counter,failed] = skipsection(fileptr,counter,target,secsize)
% Skip ahead so that counter matches target, in chunks of secsize bytes.

failed = false;
secskip = target-counter;
status = fseek(fileptr,secskip*secsize,0);
if status~=0 || feof(fileptr)
    failed = true;
end
counter = target;


