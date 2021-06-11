function CartesianRecon(pfilePath)
%% CartesianRecon - Reconstruct 2D or 3D Cartesian K-Space
%
% Copyright 2017 General Electric Company. All rights reserved.
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
    zEncoded = pfile.isZEncoded;
    clear pfile;
    
    if zEncoded == 1
        Cartesian3DRecon(pfilePath);
    else
        Cartesian2DRecon(pfilePath);
    end
