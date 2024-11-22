function result = epi_opt_param_TB(fieldmaps, rois, template, main_orientation, fov, base_res, pe_neff, delta_z, echo_spacing, TC, vx_epi, tilt, PP, dur, rfs)

% =========================================================================
% Copyright (C) 2015-2018 Steffen Volz
% Wellcome Trust Centre for Neuroimaging, London
% and Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig 
% =========================================================================

% ========================================================================= 
% This function iterates through the simulation parameter space to identify 
% the optimal BOLD (Blood Oxygen Level Dependent) sensitivity. 
% ========================================================================= 
% fieldmaps                       : Cell Array containing Field gradients 
%                                   in x,y,z directions           
% rois                            : Cell Array containing the ROIs
% template                        : Cell Array containing the Brain Mask       
% main_orientation                : Default Slice Orientation in EPI Acquistion                            
% fov                             : Field of View (in mm)
% base_res                        : Basic Resolution in pixels (Matrix size)
% pe_neff                         : Phase Encoding Steps
% delta_z                         : Slice Thickness or more precisely the
%                                   FWHM of the slice excitation profile
%                                   for a Guassian RF Pulse (in mm)
% echo_spacing                    : Echo spacing (in ms)
% TC                              : Central Echo Time (in ms)
% vx_epi                          : Voxel size (in mm)
% tilt                            : Slice angulation (in degrees)        
%                                  (1x4) array [min ref max step-size]
% PP                              : Shim gradient in z-direction (in mT/m)          
%                                  (1x4) array [min ref max step-size]  
% dur                             : Duration of compensation gradients (in ms)
%                                   (1x3) array for x,y,z directions
% rfs                             : Reduced Field Size                     
%                                   0 = no (original size), 1 = yes (1/3) 
% =========================================================================
% Updated 28/09/2024
% by Shokoufeh Golshani

% -------------------------------------------------------------------------
% Ensure PP = 0, tilt = 0 case is always evaluated
% -------------------------------------------------------------------------
if isempty(tilt)
    tilt_range = 0;
else
    tilt_range = tilt(1):tilt(4):tilt(3);
end
tilt_ref = tilt(2);

if isempty(PP)
    PP_range = 0;
else
    PP_range = PP(1):PP(4):PP(3);
end
PP_ref = PP(2);

% =========================================================================
% Phase Encoding direction 
% (0 = PE blips along y direction, 1 = PE blips along -y direction)
% =========================================================================
PE_range = 0:1;

spm_progress_bar('Init',10,'preparing ...','steps');
spm_progress_bar('Set', 0);

for n = 1:length(fieldmaps)
    display(sprintf('Using fieldmap(s) = %s;', fieldmaps{n}));
end

if (length(fieldmaps) == 3)
   fprintf('loading gradientmaps ...\n');
   rescale_gradient = 1.0*10^(-6);                           
   vol_fm_dX = spm_vol(fieldmaps{1});
   vol_fm_dY = spm_vol(fieldmaps{2});
   vol_fm_dZ = spm_vol(fieldmaps{3});
   
   fm_dX = resize(rescale_gradient*spm_read_vols(vol_fm_dX), rfs);
   fm_dY = resize(rescale_gradient*spm_read_vols(vol_fm_dY), rfs);
   fm_dZ = resize(rescale_gradient*spm_read_vols(vol_fm_dZ), rfs);

elseif (length(fieldmaps) == 1)
   fprintf('loading fieldmaps and calculating gradientmaps ...\n');   
   spm_progress_bar('Set', 1);
   [fm_dX, fm_dY, fm_dZ] = CalculateGradientmaps_TB(fieldmaps{1});
   
   fprintf('resizing gradientmaps ...\n');   
   fm_dX = resize(fm_dX, rfs);
   fm_dY = resize(fm_dY, rfs);
   fm_dZ = resize(fm_dZ, rfs);
else
  fprintf('Error: invalid number of fieldmap files;\n');
end

% -------------------------------------------------------------------------
% Setting Default Parameters
% -------------------------------------------------------------------------
fprintf('setting default scanner and EPI parameters ...\n');
spm_progress_bar('Set', 7);

epi_param_fix.main_orientation = main_orientation;
epi_param_fix.echo_spacing = echo_spacing*10^-3;

epi_param_fix.fov      = fov*10^-3;
epi_param_fix.base_res = base_res;
epi_param_fix.delta_z  = delta_z*10^-3;
epi_param_fix.TC       = TC*10^-3;
epi_param_fix.vx_epi   = vx_epi*10^-3;
epi_param_fix.TA       = echo_spacing*10^-3 * pe_neff;    % Total acquistion time

scanner_param = SetDefaultScannerParam;

% -------------------------------------------------------------------------
% Reading ROIs
% -------------------------------------------------------------------------
fprintf('reading ROIs ...\n');
spm_progress_bar('Set', 8);

for n = 1:length(rois)
    MyRoi = spm_vol(rois{n});
    ROI_slct(:,:,:,n) = resize(spm_read_vols(MyRoi), rfs);
end

% -------------------------------------------------------------------------
% Reading the template brain mask
% -------------------------------------------------------------------------
fprintf('reading the template brain mask ...\n');
spm_progress_bar('Set', 9);

for n = 1:length(template)
    tmpl_msk = spm_vol(template{n});
    Brainmask_tmpl(:,:,:,n) = resize(spm_read_vols(tmpl_msk), rfs);
end

% -------------------------------------------------------------------------
% Initializing Optimal EPI parameters
% -------------------------------------------------------------------------
PE_val = 0;
tilt_val = 0;
PP_val = 0;

epi_param_opt.GP = [0 0 PP_val]*10^-3;              % Compensation gradients (T/m)
epi_param_opt.tau = [dur(1) dur(2) dur(3)]*10^-3;   % Duration of compensation gradients (s)
epi_param_opt.tilt = tilt_val;                      % Tilt of slice (degrees)
epi_param_opt.PE_dir = PE_val;                      % PE direction (0 = standard, 1 = reversed, Siemens convention)


% -------------------------------------------------------------------------
% BS for Baseline (no tilt, no compensation, default PE direction)
% -------------------------------------------------------------------------
BS_baseline = CalculateBS_TB(fm_dX, fm_dY, fm_dZ, epi_param_opt, epi_param_fix, scanner_param);


% -------------------------------------------------------------------------
% BS for user-defined reference point
% -------------------------------------------------------------------------
if (PP_ref ~= 0) && (tilt_ref ~= 0)

    epi_param_opt.GP = [0 0 PP_ref]*10^-3;              
    epi_param_opt.tau = [dur(1) dur(2) dur(3)]*10^-3;   
    epi_param_opt.tilt = tilt_ref;                      
    epi_param_opt.PE_dir = PE_val;                      

    BS_REF = CalculateBS_TB(fm_dX, fm_dY, fm_dZ, epi_param_opt, epi_param_fix, scanner_param);

else

    BS_REF = BS_baseline;
end
% -------------------------------------------------------------------------

fprintf('exploring parameter space ... \n');

ct = 0;
ct0 = 0;

result.BS_matrix = zeros(2, size(tilt_range, 2), size(PP_range, 2), length(rois));

size_parameterspace = 2*size(tilt_range, 2)*size(PP_range, 2);

spm_progress_bar('Clear');
spm_progress_bar('Init', size_parameterspace, 'exploring parameter space', 'settings completed');
 
for PE_val = PE_range
    indx2 = 0;
    for tilt_val = tilt_range
        indx1 = 0;
        indx2 = indx2 + 1;
        for PP_val = PP_range
            indx1 = indx1 + 1;

            spm_progress_bar('Set', ct);
            ct = ct + 1;
            
            PE_all(ct) = PE_val;
            tilt_all(ct) = tilt_val;
            PP_all(ct) = PP_val;
            
            if (PE_val == 0) && (tilt_val == 0) && (PP_val == 0)
                ct0 = ct;     
            end
            
            epi_param_opt.GP = [0 0 PP_val]*10^-3;
            epi_param_opt.tau = [dur(1) dur(2) dur(3)]*10^-3;
            epi_param_opt.tilt = tilt_val;
            epi_param_opt.PE_dir = PE_val;

            BS = CalculateBS_TB(fm_dX, fm_dY, fm_dZ, epi_param_opt, epi_param_fix, scanner_param);
            
            BS_gain = (BS./(BS_baseline + eps) - 1)*100;
            
            % -------------------------------------------------------------
            % Excluding out of range values
            % -------------------------------------------------------------
            BSGainMask = (BS_gain > -100) & (BS_gain < 200);

            BrainAndGainMask = (Brainmask_tmpl > 0.99) .* BSGainMask;

            % -------------------------------------------------------------
            % Evaluating ROIS
            % -------------------------------------------------------------
            for n = 1:length(rois)

                Ind_Opt = (BrainAndGainMask.*squeeze(ROI_slct(:,:,:,n))) > 0;
                Roi_Sel_val = BS(Ind_Opt);
                Roi_mean(n, ct) = mean(Roi_Sel_val(:));
                Roi_std(n, ct) = std(Roi_Sel_val(:));
                result.BS_matrix(1 + PE_val, indx2, indx1, n) = Roi_mean(n, ct);

            end
        end
    end
end

spm_progress_bar('Clear');

fprintf('------------------------------------------------------------------------\n');
fprintf('BS Optimization done for\n');

for n = 1:length(rois)
    display(sprintf('ROI Nr. %2d: %s;', n, rois{n}));
end

fprintf(' \n');
fprintf('Optimal parameters:\n');

for n = 1:length(rois)
    [v, I] = max(squeeze(Roi_mean(n,:)));
    Ind_REF = (BS_REF.*squeeze(ROI_slct(:,:,:,n))) > 0;
    Roi_REF = mean(BS_REF(Ind_REF));
    fprintf('ROI Nr.: %2d; BS: %6.3f; BS-gain: %6.3f; BS-REF: %6.3f; PE: %1d; PP: %4.1f; tilt: %4d;\n', ...
             n, v, ((v/Roi_mean(n, ct0))-1)*100, Roi_REF, PE_all(I), PP_all(I), tilt_all(I));
    result.results(n, 1) = I;                               result.results(n, 2) = v; 
    result.results(n, 3) =  ((v/Roi_mean(n, ct0))-1)*100;   result.results(n, 4) =  PE_all(I); 
    result.results(n, 5) = PP_all(I);                       result.results(n, 6) = tilt_all(I);
    result.results(n, 7) = Roi_REF;
end

fprintf(' \n');
fprintf('BS = BS in ROI\n');
fprintf('BS-gain = BS in optimal protocol compared to default protocol\n');
fprintf('BS-REF = BS in the reference protocol\n');
fprintf('PE = phase encoding direction (0 = standard, 1 = reversed, Siemens convention)\n');
fprintf('PP = Shim gradient in z-direction\n');
fprintf('tilt = tilt of slice\n');
fprintf('------------------------------------------------------------------------\n');

% Analyse_B_matrix(result, PE_range, tilt_range, PP_range, rois)

result = 1;

end
