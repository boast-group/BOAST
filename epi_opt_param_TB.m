function result = epi_opt_param_TB(work_dir, fieldmaps, rois, main_orientation, fov, base_res, pe_neff, d, echo_spacing, TC, vx_epi, tilt_range, PP_range, res, pref)
% ========================================================================
% loop through protocol parameter space and optimize BOLD sensitivity
% Copyright (C) 2015-2018 Steffen Volz
% Wellcome Trust Centre for Neuroimaging, London
% and Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig 
% ========================================================================

%reducing size of maps(1mm) by factor of 3
%res = 1;

spm_progress_bar('Init',10,'preparing ...','steps');
spm_progress_bar('Set',0);
for n = 1:length(fieldmaps)
display(sprintf('Using fieldmap(s) = %s;',fieldmaps{n}));
end

if (length(fieldmaps) == 3)
   display(sprintf('loading gradientmaps ...'));
   rescale_gradient = 1.0;
   vol_fm_dX = spm_vol(fieldmaps{1});
   vol_fm_dY = spm_vol(fieldmaps{2});
   vol_fm_dZ = spm_vol(fieldmaps{3});
   fm_dX = myresize_3(rescale_gradient*0.000001*spm_read_vols(vol_fm_dX),res);
   fm_dY = myresize_3(rescale_gradient*0.000001*spm_read_vols(vol_fm_dY),res);
   fm_dZ = myresize_3(rescale_gradient*0.000001*spm_read_vols(vol_fm_dZ),res);
elseif (length(fieldmaps) == 1)
   display(sprintf('loading fieldmaps and calculating gradientmaps ...'));
   spm_progress_bar('Set',1);
   [fm_dX,fm_dY,fm_dZ] = CalculateGradientmaps_TB(fieldmaps{1});
   display(sprintf('resizing gradientmaps ...'));
   fm_dX = myresize_3(fm_dX,res);
   fm_dY = myresize_3(fm_dY,res);
   fm_dZ = myresize_3(fm_dZ,res);
else
%   display(sprintf('Error: invalid number of fieldmap files;'));
end

display(sprintf('loading default scanner and EPI parameter ...'));
spm_progress_bar('Set',7);
% load default scanner and EPI parameter
scanner_param = SetDefaultScannerParam;
epi_param_fix = SetDefaultEPIParam;

display(sprintf('setting user defined parameter ...'));
% set user defined parameter
epi_param_fix.main_orientation = main_orientation;
epi_param_fix.echo_spacing = echo_spacing; % echo spacing in ms

epi_param_fix.fov      = fov;
epi_param_fix.base_res = base_res;
epi_param_fix.d        = d;
epi_param_fix.TC       = TC;
epi_param_fix.vx_epi   = vx_epi;
epi_param_fix.TA       = epi_param_fix.echo_spacing*pe_neff;

display(sprintf('reading ROIs ...'));
spm_progress_bar('Set',8);
for n = 1:length(rois)
%   display(sprintf('Reading ROI: %s;',rois{n}));
   MyRoi = spm_vol(rois{n});
   Roi_Sel(:,:,:,n) = myresize_3(spm_read_vols(MyRoi),res);
end
spm_progress_bar('Set',9);

% this solution is not extremely elegant
Brainmask = squeeze(Roi_Sel(:,:,:,7));

% store BS for default sequence
PE_val=0;tilt_val=0;PP_val=0;
%epi_param_opt.TC = TC; % echo time in protocol
epi_param_opt.GP = [0 0 PP_val]*10^-6; % compensation gradient
epi_param_opt.a = tilt_val;%  a: tilt of slice (deg, T>C, Siemens convention)
epi_param_opt.PE_dir = PE_val;% PE_dir: PE direction of EPI (0 = standard, 1 = reversed, Siemens convention)
BS0 = CalculateBS_TB(0, fm_dX, fm_dY, fm_dZ, epi_param_opt, epi_param_fix, scanner_param);

display(sprintf('exploring parameter space ... '));

ct = 0;
ct0 = 0;

result.BS_matrix = zeros(2,size(tilt_range,2), size(PP_range,2),7);

size_parameterspace = 2*size(tilt_range,2)*size(PP_range,2);
spm_progress_bar('Clear');
spm_progress_bar('Init',size_parameterspace,'exploring parameter space','settings completed');

for PE_val = 0:1:1
 for tilt_val = tilt_range
  for PP_val = PP_range
   spm_progress_bar('Set',ct);
   ct = ct+1;
 
   PE_all(ct)=PE_val;
   tilt_all(ct)=tilt_val;
   PP_all(ct)=PP_val;
  
   if (PE_val==0) & (tilt_val==0) & (PP_val==0)
      ct0 = ct;
   end
  
   %epi_param_opt.TC = TC; % echo time in protocol
   epi_param_opt.GP = [0 0 PP_val]*10^-6; % compensation gradient
   epi_param_opt.a = tilt_val;%  a: tilt of slice (deg, T>C, Siemens convention)
   epi_param_opt.PE_dir = PE_val;% PE_dir: PE direction of EPI (0 = standard, 1 = reversed, Siemens convention)

   BS = CalculateBS_TB(0, fm_dX, fm_dY, fm_dZ, epi_param_opt, epi_param_fix, scanner_param);

   BS_gain = (BS./(BS0+eps)-1)*100;
   BSGainMask = (BS_gain > -100) & (BS_gain < 200);% for excluding stupid values
   BrainAndGainMask = (Brainmask > 0.99).*BSGainMask;
  
% evaluate ROIS
   for n = 1:length(rois)
     Ind_ToOpt = find((BrainAndGainMask.*squeeze(Roi_Sel(:,:,:,n))) > 0);
     Roi_Sel_val = BS(Ind_ToOpt);
     Roi_mean(n,ct) = mean(Roi_Sel_val(:));
     Roi_std(n,ct)  = std(Roi_Sel_val(:));
     result.BS_matrix(1+PE_val, 1+(tilt_val-tilt_range(1))/((tilt_range(size(tilt_range,2))-tilt_range(1))/(size(tilt_range,2)-1)), 1+(PP_val-PP_range(1))/((PP_range(size(PP_range,2))-PP_range(1))/(size(PP_range,2)-1)),n) = Roi_mean(n, ct);
   end
  
  end
 end
end

spm_progress_bar('Clear');

display(sprintf('------------------------------------------------------------------------'));
display(sprintf('BS Optimization done for'));

for n = 1:length(rois)
display(sprintf('ROI Nr. %2d: %s;',n,rois{n}));
end

display(sprintf(' '));
display(sprintf('Optimal parameters:'));
for n = 1:length(rois)
   [v I] = max(squeeze(Roi_mean(n,:)));
   display(sprintf('ROI Nr.: %2d; BS: %6.3f; BS-gain: %6.3f; PE: %1d; PP: %4.1f; tilt: %4d;',n, v, v/Roi_mean(n,ct0), PE_all(I), PP_all(I), tilt_all(I)));
   result.results(n,1) = I; result.results(n,2) = v; result.results(n,3) =  v/Roi_mean(n,ct0); result.results(n,4) =  PE_all(I); result.results(n,5) = PP_all(I); result.results(n,6) = tilt_all(I);
end

display(sprintf(' '));
display(sprintf('BS = BS in ROI compared to BS zero-gradients'));
display(sprintf('BS-gain = BS in optimal protocol compared to default protocol'));
display(sprintf('PE = phase encoding direction (0 = standard, 1 = reversed, Siemens convention)'));
display(sprintf('PP = Shim gradient in z-direction'));
display(sprintf('tilt = tilt of slice'));
display(sprintf('------------------------------------------------------------------------'));

result = 1;

end
