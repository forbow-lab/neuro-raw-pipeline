function pushDICOMs(folder)

addpath('/biotic/opt/matlab/toolbox/orchestra');
networkHandle = GERecon('Dicom.CreateNetwork', '142.239.156.234', 4006, 'HIDIRES', 'CHDIAW01');
GERecon('Dicom.Store', folder, networkHandle);
GERecon('Dicom.CloseNetwork', networkHandle);

end

