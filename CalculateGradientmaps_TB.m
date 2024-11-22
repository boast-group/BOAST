function [fm_dX, fm_dY, fm_dZ] = CalculateGradientmaps_TB(file_fmp) 

% =========================================================================
% Copyright (C) 2014-2018 Steffen Volz
% Wellcome Trust Centre for Neuroimaging, London
% and Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig 
% =========================================================================

% ========================================================================= 
% This function calculates the field map gradients in the x,y,z directions 
% from a field map with 1mm isotropic voxels.
% ========================================================================= 
% Input:
%        file_fmp                 : Field map file obtained from the fieldmap 
%                                   toolbox (in Hz), with a resolution of 1 mm.
% Output:
%        fm_dX, fm_dY, fm_dZ      : Field Gradients in x,y,z directions in T/m 
% =========================================================================                                   in x,y,z directions           
% Updated 29/09/2024
% By Shokoufeh Golshani


gam = 2*pi*42.58*10^6;        % gyromagnetic ratio in Hz/Tesla for protons
Larmor = gam/2/pi;

% =========================================================================
% Read field map and calculate the gradient of the field
% =========================================================================
fm_V = spm_vol(file_fmp);
rot_mat = fm_V(1).mat(1:3,1:3)*10^-3;
vx_size = [1 1 1].*10^(-3);    

if diag(vx_size(1))~=abs(rot_mat(1))
    error('CalculateGradientmaps: field map seems to be spatially transformed!');
end

spm_progress_bar('Set', 1);

[x, y, z] = ndgrid(1:fm_V.dim(1), 1:fm_V.dim(2), 1:fm_V.dim(3));

spm_progress_bar('Set', 3);

[fm_val, fm_dX, fm_dY, fm_dZ] = spm_sample_vol(fm_V, x, y, z, -8);
fm_val = reshape(fm_val,[fm_V.dim(1), fm_V.dim(2), fm_V.dim(3)]);

fm_dX = reshape(fm_dX, [fm_V.dim(1), fm_V.dim(2), fm_V.dim(3)]);
fm_dY = reshape(fm_dY, [fm_V.dim(1), fm_V.dim(2), fm_V.dim(3)]);
fm_dZ = reshape(fm_dZ, [fm_V.dim(1), fm_V.dim(2), fm_V.dim(3)]);

spm_progress_bar('Set', 6);
fm_dX = fm_dX/vx_size(1)/Larmor;
fm_dY = fm_dY/vx_size(2)/Larmor;
fm_dZ = fm_dZ/vx_size(3)/Larmor;

end
