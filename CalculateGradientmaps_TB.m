function [fm_dX, fm_dY, fm_dZ] = CalculateGradientmaps_TB(file_fmp) 

% =========================================================================
% Calculate the fieldmap gradients in three x,y,z directions
% Copyright (C) 2014-2018 Steffen Volz
% Wellcome Trust Centre for Neuroimaging, London
% and Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig 
% =========================================================================

% Updated 29/09/2024
% By Shokoufeh Golshani

gam = 2*pi*42.58*10^6;        % gyromagnetic ratio in Hz/Tesla for protons
Larmor = gam/2/pi;

% =========================================================================
% Read field map and calculate the gradient of the field
% =========================================================================
fm_V = spm_vol(file_fmp);
rot_mat = fm_V(1).mat(1:3,1:3)*10^-3;
vx_size = [0.001, 0.001, 0.001];

if diag(vx_size)~=abs(rot_mat)
    error('CalculateGradientmaps: field map seems to be spatially transformed!');
end

spm_progress_bar('Set', 1);

[x,y,z] = ndgrid(1:fm_V.dim(1), 1:fm_V.dim(2), 1:fm_V.dim(3));

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
