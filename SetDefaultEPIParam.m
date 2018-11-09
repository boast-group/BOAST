function epi_param = SetDefaultEPIParam

our_epi.main_orientation = 'TRA';
our_epi.fov = 0.192; % field of view in EPI 192 mm
our_epi.base_res = 64;
% slice thickness
our_epi.d=2*10^-3; %Siemens pulse approximates Gaussian with 2 mm FWHM
%our_epi.echo_spacing = 0.330; % for allegra
our_epi.echo_spacing = 0.500; % echo spacing in ms
our_epi.TC= 30; % ms
our_epi.vx_epi=[3 3 3];

our_epi.vx_epi = our_epi.vx_epi*10^-3;
our_epi.echo_spacing = our_epi.echo_spacing*10^-3;
our_epi.TC = our_epi.TC*10^-3;
our_epi.TA = our_epi.echo_spacing*our_epi.base_res;

epi_param = our_epi;

end
