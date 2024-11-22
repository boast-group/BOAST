function BS = CalculateBS_TB(fm_dX, fm_dY, fm_dZ, epi_param_opt, epi_param_fix, scanner_param) 

% =========================================================================
% This function calculates BOLD sensitivity using field map gradients and 
% a defined set of parameters.
% Copyright (C) 2014-2018 Steffen Volz
% Wellcome Trust Centre for Neuroimaging, London
% and Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig 
% =========================================================================

% Updated 28/09/2024
% By Shokoufeh Golshani

% =========================================================================
% Default variables
% =========================================================================
gam = 2*pi*42.58*10^6;         % gyromagnetic ratio in Hz/Tesla for protons

  
% =========================================================================
% Compensation Gradients
% =========================================================================
GPrep_RO = epi_param_opt.GP(1);             tau_RO = epi_param_opt.tau(1);
GPrep_PE = epi_param_opt.GP(2);             tau_PE = epi_param_opt.tau(2);
GPrep_S = epi_param_opt.GP(3);              tau_S = epi_param_opt.tau(3);

Angle = -epi_param_opt.tilt/180*pi;

% =========================================================================
% Rotation matrix for converting gradient from XYZ to PRS
% Refer to the Fig. 1 in the paper for more explanation (TRA > COR, Siemens convention)
% =========================================================================
if strcmp(epi_param_fix.main_orientation,'TRA') == 1
    pe_vec = [0 cos(Angle) sin(Angle)];
    ro_vec = [1 0 0];
    sl_vec = [0 sin(Angle) -cos(Angle)];
elseif strcmp(epi_param_fix.main_orientation,'SAG') == 1
    pe_vec = [0 cos(Angle) sin(Angle)];
    ro_vec = [0 -sin(Angle) cos(Angle)];
    sl_vec = [1 0 0];    
elseif strcmp(epi_param_fix.main_orientation,'COR') == 1
    pe_vec = [0 -sin(Angle) cos(Angle)];
    ro_vec = [1 0 0];
    sl_vec = [0 cos(Angle) sin(Angle)];
end
    
if epi_param_opt.PE_dir == 1                        % k-space down
    pe_vec= -pe_vec;
end

fGP = fm_dX*pe_vec(1) + fm_dY*pe_vec(2) + fm_dZ*pe_vec(3);
fGR = fm_dX*ro_vec(1) + fm_dY*ro_vec(2) + fm_dZ*ro_vec(3);
fGS = fm_dX*sl_vec(1) + fm_dY*sl_vec(2) + fm_dZ*sl_vec(3);

% =========================================================================
% Calculate the Q value which determines distortion and echo shift
% =========================================================================
Q = 1 - (gam * epi_param_fix.echo_spacing/2/pi * epi_param_fix.fov * fGP);


% =========================================================================
% Actual Echo time
% =========================================================================
TE_shift = gam * GPrep_PE * tau_PE/2/pi * epi_param_fix.fov * epi_param_fix.echo_spacing;

TE = (epi_param_fix.TC + TE_shift)./Q;
dTE = (TE - epi_param_fix.TC); 

% =========================================================================
% Sudden dropout when echo is shifted out of acq window in the PE direction
% =========================================================================
shift_mask = ones(size(Q));

shift_mask(abs(dTE) > (epi_param_fix.TA/2)) = 0;
shift_mask = shift_mask.*(abs(fGR + GPrep_RO.*tau_RO./TE) < (pi/gam./TE/epi_param_fix.vx_epi(1)));

% =========================================================================
% Through-plane gradient --- Gaussian RF pulse
% =========================================================================
I = exp((-gam^2/16/log(2)*epi_param_fix.delta_z^2).*((GPrep_S.*tau_S + fGS.*TE).^2));

% =========================================================================
% BOLD Sensitivity
% =========================================================================
BS = I./Q.^2.*exp(-epi_param_fix.TC/scanner_param.T2s.*(1./Q-1));

% =========================================================================
% Correct maps for shifts of data out of acquisition window
% =========================================================================
BS = BS.*shift_mask;


end
