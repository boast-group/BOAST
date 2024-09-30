function scanner_param = SetDefaultScannerParam

% Updated 23/09/2024
% by Shokoufeh Golshani

% =========================================================================
% This function sets the scanner parameters. 
% These values can be modified as needed to suit specific requirements.
% =========================================================================
% name                            : Scanner name
% B0                              : Field strength
% T2s                             : T2* value - An averaged value over the 
%                                   entire brain or a voxel-wise map
%                                   Although not directly a scanner parameter, 
%                                   it depends on the field strength!
% =========================================================================

fixedScannerparam.name = 'Trio';
fixedScannerparam.B0 = 3;

answer = questdlg('Do you have a T2* map?', ...
	'T2* impact', ...
	'Yes','No', 'No');

switch answer
    case 'Yes'
        [file, location] = uigetfile({'*.nii'}, 'File Selector');
        fixedScannerparam.T2s = double(niftiread(fullfile(location, file)));
    
    case 'No'
        % T2* at 3 T (Wansapura et al., JMRI 1999)
        fixedScannerparam.T2s = 45*10^-3;  
end

scanner_param = fixedScannerparam;

end
