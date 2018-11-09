function BS = CalculateBS_TB(default, fm_dX, fm_dY, fm_dZ, epi_param_opt, epi_param_fix, scanner_param) 
% ========================================================================
% calculate BOLD sensitivity from gradient fieldmap
% based on calc_BS_fm_atlas (NW)
% Copyright (C) 2014-2018 Steffen Volz
% Wellcome Trust Centre for Neuroimaging, London
% and Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig 
% ========================================================================

if default==1
    scanner_param = SetDefaultScannerParam;
    epi_param_fix = SetDefaultEPIParam;
end
% else take epi_param_fix, scanner_param as is

gam = 2*pi*42580000; % gyromagnetic ratio 2*pi*42.58 MHz/Tesla for protons
Larmor = gam/2/pi;

GPS = epi_param_opt.GP(3);
GPRO = epi_param_opt.GP(1);
GPPE = epi_param_opt.GP(2);

a=-epi_param_opt.a/180*pi;

% rotation matrix for converting gradient from XYZ to PRS
if strcmp(epi_param_fix.main_orientation,'TRA')==1
    pe_vec = [0 cos(a) sin(a)];
    ro_vec = [1 0 0];
    sl_vec = [0 sin(a) -cos(a)];
elseif strcmp(epi_param_fix.main_orientation,'SAG')==1
% 2. for SAG slice (RO HF, PE PA, SL RL) tilted inplane angle a
    pe_vec = [0 cos(a) sin(a)];
    ro_vec = [0 -sin(a) cos(a)];
    sl_vec = [1 0 0];    
elseif strcmp(epi_param_fix.main_orientation,'COR')==1
%{
% 3b. for COR slice (RO HF, PE RL, SL AP) tilted to COR angle a
pe_vec = [1 0 0];
ro_vec = [0 -sin(a) cos(a)];
sl_vec = [0 -cos(a) -sin(a)];
%}

% 3a. for COR slice (RO RL, PE HF, SL PA) tilted inplane angle a
    pe_vec = [0 -sin(a) cos(a)];
    ro_vec = [1 0 0];
    sl_vec = [0 cos(a) sin(a)];
end
    
if epi_param_opt.PE_dir == 1 % k space down
 pe_vec= -pe_vec;
end

fGP = fm_dX*pe_vec(1)+fm_dY*pe_vec(2)+fm_dZ*pe_vec(3);
fGR = fm_dX*ro_vec(1)+fm_dY*ro_vec(2)+fm_dZ*ro_vec(3);
fGS = fm_dX*sl_vec(1)+fm_dY*sl_vec(2)+fm_dZ*sl_vec(3);

% Calculate the Q value which determines distortion and echo shift
Q = 1-(gam*epi_param_fix.echo_spacing/2/pi*epi_param_fix.fov*fGP);
% Compensation gradient in PE direction
GPPE_shift = gam * GPPE /2/pi * epi_param_fix.fov * epi_param_fix.echo_spacing;
% real echo times
TE = (epi_param_fix.TC+GPPE_shift)./Q;
dTE = (TE-epi_param_fix.TC); % actual TE
 
shift_mask = ones(size(Q));
PE_shift_mask = ones(size(Q));
RO_shift_mask = ones(size(Q));

% sudden dropout when echo is shifted out of acq window in the PE direction
shift_mask(abs(dTE)>(epi_param_fix.TA/2))=0;
PE_shift_mask(abs(dTE)>(epi_param_fix.TA/2))=0;
% sudden dropout when echo is shifted out of acq window in the RO direction
% shift_mask = shift_mask.*(abs(fm_dX)<(pi/gam/TE/vx_epi(1)));

shift_mask = shift_mask.*(abs(fGR+GPRO./TE)<(pi/gam./TE/epi_param_fix.vx_epi(1)));
RO_shift_mask = RO_shift_mask.*(abs(fGR+GPRO./TE)<(pi/gam./TE/epi_param_fix.vx_epi(1)));

% Through-plane gradient
% --- Gaussian RF pulse
I = exp(-gam^2/16/log(2)*epi_param_fix.d^2*((GPS+fGS.*TE).^2)); % Gaussian slice profile

BS = I./Q.^2.*exp(-epi_param_fix.TC/scanner_param.T2s.*(1./Q-1));

% correct maps for shifts of data out of acquisition window
BS=BS.*shift_mask;

I_real = I./Q.*exp(-epi_param_fix.TC/scanner_param.T2s.*(1./Q-1)); % correct signal for echo shifting
I_real=I_real.*shift_mask;

result = BS;
end

