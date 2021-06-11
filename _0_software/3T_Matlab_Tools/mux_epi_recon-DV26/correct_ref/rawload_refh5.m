function pha_coe = rawload_refh5(h5file, ny_pres, nslices, ncoils, frames, slices, coils)



if ~exist('frames', 'var');      frames = [];            end;
if ~exist('slices', 'var');      slices = [];            end;
if ~exist('coils', 'var');       coils  = [];            end;

if (isempty(frames));            frames = 1:ny_pres;         end;
if (isempty(slices));            slices = 1:nslices;         end;
if (isempty(coils));             coils  = 1:ncoils;          end;



% pha_coe = zeros(2, ny_pres, nslices, ncoils);
pha_coe = zeros(2, frames, nslices, ncoils);

for sl = 1:nslices, 
    for ch = 1:ncoils, 
        pha_coe(1, :, sl, ch) = h5read(h5file, sprintf('/slice_%d/channel_%d/ConstantCoefficients', sl-1, ch-1));
        pha_coe(2, :, sl, ch) = h5read(h5file, sprintf('/slice_%d/channel_%d/LinearCoefficients', sl-1, ch-1));
    end
end

pha_coe = pha_coe(:, 1:frames, slices, coils);
