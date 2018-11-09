function fmpoptbs = tbx_cfg_fmpoptbs
% MATLABBATCH Configuration file for toolbox 'tbx_cfg_FmpOptBS'
%_______________________________________________________________________
% Copyright (C) 2015-2018 Steffen Volz
% Wellcome Trust Centre for Neuroimaging, London
% and Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig 

% $Id: tbx_cfg_FmpOptBS.m 89 2015-12-16 11:03:00Z steffen $

if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','FmpOptBS')); end

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Default values that are common to all tbx_cfg_FmpOptBS jobs
% ---------------------------------------------------------------------
%

% ---------------------------------------------------------------------
% field maps
% ---------------------------------------------------------------------
fieldmaps         = cfg_files;
fieldmaps.tag     = 'fieldmaps';
fieldmaps.name    = 'Input fieldmaps';
fieldmaps.val{1}  = {'files/mean_wr_calcDX.nii' 'files/mean_wr_calcDY.nii' 'files/mean_wr_calcDZ.nii'};
%fieldmaps.val{1}  = {'files/mean_wr_fmp.nii'};
fieldmaps.help    = {'Input one fieldmap or 3 fieldmap gradient (dX dY dZ) files for optimizing BOLD sensitivity.'};
%fieldmaps.filter  = 'fmp';
fieldmaps.ufilter = '.*';
fieldmaps.num     = [1 Inf];
% ---------------------------------------------------------------------
% template
% ---------------------------------------------------------------------
template         = cfg_files;
template.tag     = 'template';
template.name    = 'Input template';
template.val{1}  = {'/files/Brainmask.nii'};
template.help    = {'template for illustration (only placeholder at the moment)'};
template.ufilter = '.*';
template.num     = [1 Inf];
% ---------------------------------------------------------------------
% ROIs
% ---------------------------------------------------------------------
rois         = cfg_files;
rois.tag     = 'rois';
rois.name    = 'ROI files';
rois.val{1}  = {'files/MyRoi_mOFC-and-rACC_mod.nii' 'files/MyRoi_InferiorTemporalGyri.nii' 'files/MyRoi_TemporalPoles.nii' 'files/MyRoi_Amygdala.nii' 'files/MyRoi_Hippo-andParahippocampalRegion.nii' 'files/MyRoi_WellShimmed.nii' 'files/Brainmask.nii'};
rois.help    = {'ROIs for which BOLD is optimized.'};
rois.ufilter = '.*';
rois.num     = [1 Inf];
% ---------------------------------------------------------------------
% Input Files
% ---------------------------------------------------------------------
inputfiles         = cfg_branch;
inputfiles.tag     = 'inputfiles';
inputfiles.name    = 'Input Files';
inputfiles.val     = {fieldmaps template rois};
inputfiles.help    = {'Needed Input Files'};
% ---------------------------------------------------------------------
% menu main orientation
% ---------------------------------------------------------------------
main_orientation         = cfg_menu;
main_orientation.tag     = 'main_orientation';
main_orientation.name    = 'Choose the main orientation';
main_orientation.help    = {'Option to choose the main orientation'};
main_orientation.labels  = {'TRA' 'COR' 'SAG'};
main_orientation.values  = {'TRA' 'COR' 'SAG'};
main_orientation.val     = {'TRA'};
% ---------------------------------------------------------------------
% Field of View
% ---------------------------------------------------------------------
fov         = cfg_entry;
fov.tag     = 'fov';
fov.name    = 'field of view';
fov.val     = {192};
fov.help    = {'field of view in mm'};
fov.strtype = 'r';
fov.num     = [1 1];
% ---------------------------------------------------------------------
% Base Resolution
% ---------------------------------------------------------------------
base_res         = cfg_entry;
base_res.tag     = 'base_res';
base_res.name    = 'Base Resolution';
base_res.val     = {64};
base_res.help    = {'Base Resolution in #px'};
base_res.strtype = 'r';
base_res.num     = [1 1];
% ---------------------------------------------------------------------
% Phase Encoding Steps
% ---------------------------------------------------------------------
pe_neff         = cfg_entry;
pe_neff.tag     = 'pe_neff';
pe_neff.name    = 'Phase Encoding Steps';
pe_neff.val     = {72};
pe_neff.help    = {'Phase Encoding Steps in #px'};
pe_neff.strtype = 'r';
pe_neff.num     = [1 1];
% ---------------------------------------------------------------------
% Slice Thickness
% ---------------------------------------------------------------------
slicethickness         = cfg_entry;
slicethickness.tag     = 'slicethickness';
slicethickness.name    = 'Slice Thickness';
slicethickness.val     = {3};
slicethickness.help    = {'Slice Thickness in mm'};
slicethickness.strtype = 'r';
slicethickness.num     = [1 1];
% ---------------------------------------------------------------------
% Echo Spacing
% ---------------------------------------------------------------------
echospacing         = cfg_entry;
echospacing.tag     = 'echospacing';
echospacing.name    = 'Echo Spacing';
echospacing.val     = {0.5};
echospacing.help    = {'Echo Spacing in ms'};
echospacing.strtype = 'r';
echospacing.num     = [1 1];
% ---------------------------------------------------------------------
% Echo Time
% ---------------------------------------------------------------------
echotime         = cfg_entry;
echotime.tag     = 'echotime';
echotime.name    = 'Echo Time';
echotime.val     = {30};
echotime.help    = {'Echo Time in ms'};
echotime.strtype = 'r';
echotime.num     = [1 1];
% ---------------------------------------------------------------------
% vox Voxel size
% ---------------------------------------------------------------------
vox         = cfg_entry;
vox.tag     = 'vox';
vox.name    = 'Voxel size';
vox.val     = {[3 3 3]};
vox.help    = {'voxel size [phase, read, slice] in mm'};
vox.strtype = 'r';
vox.num     = [1 3];
% ---------------------------------------------------------------------
% Fixed Protocol Parameters
% ---------------------------------------------------------------------
fixedparameters         = cfg_branch;
fixedparameters.tag     = 'fixedparameters';
fixedparameters.name    = 'Fixed Protocol Parameters';
fixedparameters.val     = {main_orientation fov base_res pe_neff slicethickness echospacing echotime vox};
fixedparameters.help    = {'Fixed Protocol Parameters'};
% ---------------------------------------------------------------------
% Parameter shimz
% ---------------------------------------------------------------------
shimz         = cfg_entry;
shimz.tag     = 'shimz';
shimz.name    = 'shimz';
shimz.val     = {[-3 0 3 0.5]};
shimz.help    = {'Shim gradient in z-direction [min ref max step-size]'};
shimz.strtype = 'r';
shimz.num     = [1 4];
% ---------------------------------------------------------------------
% Parameter tilt
% ---------------------------------------------------------------------
tilt         = cfg_entry;
tilt.tag     = 'tilt';
tilt.name    = 'tilt';
tilt.val     = {[-45 0 45 5]};
tilt.help    = {'tilt in degree [min ref max step-size]'};
tilt.strtype = 'r';
tilt.num     = [1 4];
% ---------------------------------------------------------------------
% Simulation Parameters
% ---------------------------------------------------------------------
simu         = cfg_branch;
simu.tag     = 'simu';
simu.name    = 'Simulation Parameters';
simu.val     = {shimz tilt};
%simu.val     = {TE,PEdir,tilt,shimx,shimy,shimz};
simu.help    = {'Parameters to be simulated: all these parameters have a minimum, maximum and default value and a step size for the optimization procedure'};
% ---------------------------------------------------------------------
% Reduce Field size
% ---------------------------------------------------------------------
rfs         = cfg_entry;
rfs.tag     = 'rfs';
rfs.name    = 'Reduce Field size';
rfs.val     = {1};
rfs.help    = {'0 = no (original size), 1 = yes (1/3)'};
rfs.strtype = 'r';
rfs.num     = [1 1];
% ---------------------------------------------------------------------
% Other Settings
% ---------------------------------------------------------------------
other         = cfg_branch;
other.tag     = 'other';
other.name    = 'Other Settings';
other.val     = {rfs};
other.help    = {'Other Settings used for optimization'};
% ---------------------------------------------------------------------
% preproc8 Segment
% ---------------------------------------------------------------------
fmpoptbs         = cfg_exbranch;
fmpoptbs.tag     = 'FmpOptBS';
fmpoptbs.name    = 'BS Optimisation';
fmpoptbs.val     = {inputfiles fixedparameters simu other};
fmpoptbs.help    = {'This toolbox is currently only work in progress.'};
fmpoptbs.prog = @fmpoptbs_apply;
fmpoptbs.vout = @vout_fmpoptbs_apply;
%----------------------------------------------------------------------

function out = fmpoptbs_apply(job)

%length(job.inputfiles.fieldmaps)
%for n = 1:length(job.inputfiles.fieldmaps)
%display(sprintf('fieldmaps = %s;',job.inputfiles.fieldmaps{n}));
%end
%length(job.inputfiles.rois)
%for n = 1:length(job.inputfiles.rois)
%display(sprintf('rois = %s;',job.inputfiles.rois{n}));
%end

%display(sprintf('fixed parameters: main_orientation = %s, fov = %5.1f, base_res = %5.1f, slicethickness = %5.1f, echospacing = %5.3f, echotime = %5.1f;', job.fixedparameters.main_orientation, job.fixedparameters.fov, job.fixedparameters.base_res, job.fixedparameters.slicethickness, job.fixedparameters.echospacing, job.fixedparameters.echotime));%vox
%display(sprintf('shimz from %5.1f to %5.1f in steps of %5.1f with reference %5.1f;',job.simu.shimz(1),job.simu.shimz(3),job.simu.shimz(4),job.simu.shimz(2)));
%display(sprintf('tilt from %5.1f to %5.1f in steps of %5.1f with reference %5.1f;',job.simu.tilt(1),job.simu.tilt(3),job.simu.tilt(4),job.simu.tilt(2)));

TB_files_dir = sprintf('%s/FmpOptBS/files/', spm('dir'));
opt_results = epi_opt_param_TB(TB_files_dir, job.inputfiles.fieldmaps, job.inputfiles.rois, job.fixedparameters.main_orientation, job.fixedparameters.fov*10^-3, job.fixedparameters.base_res, job.fixedparameters.pe_neff, job.fixedparameters.slicethickness*10^-3, job.fixedparameters.echospacing*10^-3, job.fixedparameters.echotime*10^-3, job.fixedparameters.vox*10^-3, job.simu.tilt(1):job.simu.tilt(4):job.simu.tilt(3), job.simu.shimz(1):job.simu.shimz(4):job.simu.shimz(3), job.other.rfs, 'Opt_');


out.fmfiles = job.inputfiles.fieldmaps;
function dep = vout_fmpoptbs_apply(job)
% do something
dep = cfg_dep;

% end;
%------------------------------------------------------------------------
