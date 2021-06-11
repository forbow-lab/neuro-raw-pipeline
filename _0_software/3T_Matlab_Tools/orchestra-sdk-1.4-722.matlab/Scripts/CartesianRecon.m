function CartesianRecon(pfilePath)
%% CartesianRecon - Reconstruct 2D Cartesian K-Space
%
% Copyright 2014 General Electric Company. All rights reserved.
% GE Proprietary and Confidential Information. Only to be distributed with
% permission from GE. Resulting outputs are not for diagnostic purposes.
%
% CartesianRecon(pfilePath)
% will reconstruct the Cartesian K-Space in the given pfile. Except for
% 3D scans with ARC, this includes Pfiles from both 2D and 3D
% acquisitions.
%
% Limitations: Parallel imaging, intensity correction

    % Load Pfile
    pfile = GERecon('Pfile.Load', pfilePath);

    for p = 1:pfile.phases
        for s = 1:pfile.slices

            % Get corners and orientation for this slice location
            corners = GERecon('Pfile.Corners', s);
            orientation = GERecon('Pfile.Orientation', s);

            for e = 1:pfile.echoes
                for c = 1:pfile.channels

                    % Load K-Space
                    kspace = GERecon('Pfile.KSpace', s, e, c);

                    % Transform K-Space
                    channelImage = GERecon('Transform', kspace);

                    channelImages(:,:,c) = channelImage;
                end

                % Apply Channel Combination
                combinedImage = GERecon('SumOfSquares', channelImages);

                % Create Magnitude Image
                magnitudeImage = abs(combinedImage);

                % Apply Gradwarp
                gradwarpedImage = GERecon('Gradwarp', magnitudeImage, corners);

                % Orient the image
                finalImage = GERecon('Orient', gradwarpedImage, orientation);

                % Display
                imagesc(finalImage);
                title(['Phase: ' num2str(p) 'Slice: ' num2str(s) 'Echo: ' num2str(e)]);

                % Save DICOMs
                imageNumber = ImageNumber(s, e, p, pfile);
                filename = ['DICOMs/image' num2str(imageNumber) '.dcm'];
                GERecon('Dicom.Write', filename, finalImage, imageNumber, orientation, corners);

                pause(0.1);
            end
        end
    end
end

function number = ImageNumber(slice, echo, phase, pfile)
% Image numbering scheme:
% P0S0E0, P0S0E1, ... P0S0En, P0S1E0, P0S1E1, ... P0S1En, ... P0SnEn, ...
% P1S0E0, P1S0E1, ... PnSnEn
    slicesPerPhase = pfile.slices * pfile.echoes;
    number = (phase-1) * slicesPerPhase + (slice-1) * pfile.echoes + (echo-1) + 1;
end
