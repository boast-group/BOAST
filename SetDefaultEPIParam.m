function epi_param = SetDefaultEPIParam

% Updated 23/09/2024
% by Shokoufeh Golshani

% =========================================================================
% This function sets the fixed parameters for the EPI protocol. 
% These values can be modified as needed to suit specific requirements 
% or preferences.
% =========================================================================
% main_orientation                : Slice oriantation
%                                   'TRA' : transverse 
%                                   'CRO' : coronal
%                                   'SAG' : sagittal
% fov                             : Field of view (in mm)
% base_res                        : Basic resolution (Matrix size)
% pe_neff                         : Effective phase encoding steps in
%                                   Siemense scanner
% delta_z                         : Full width at half-maximum (FWHM) 
%                                   of the slice profile (in mm)
% echo_spacing                    : Echo spacing (in ms)
% echotime                        : Effective (central) echo time (in ms)
% vox                             : voxel size (in mm) 
%                                   1x3 array (x y z direction)
% =========================================================================

fixedEPIparam.main_orientation = 'TRA';
fixedEPIparam.fov = 192;                         
fixedEPIparam.base_res = 64;
fixedEPIparam.pe_neff = 72;
fixedEPIparam.delta_z = 2;         % FWHM for a Gaussian profile in Siemens 
                                   % scanner using slice thickness = 3 mm
fixedEPIparam.echo_spacing = 0.5;
fixedEPIparam.echotime = 30; 
fixedEPIparam.vox = [3 3 3];

epi_param = fixedEPIparam;

end
