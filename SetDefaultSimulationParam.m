function sim_param = SetDefaultSimulationParam

% Created 23/09/2024
% by Shokoufeh Golshani

% ========================================================================= 
% This function sets the desired simulation parameters. 
% These values can be modified as needed to suit specific requirements.
% ========================================================================= 
% shimz                           : Shim gradient in z-direction (in mT/m)          
%                                  (1x4) array [min ref max step-size]  
% tau                             : Duration of compensation gradients (in ms)
%                                   (1x3) array for x,y,z directions
% tilt                            : Slice angulation (in degrees)        
%                                  (1x4) array [min ref max step-size]
% rfs                             : Reduced Field Size                     
%                                   0 = no (original size), 1 = yes (1/3) 
% =========================================================================

SimuParam.shimz = [-5 0 5 0.5];
SimuParam.tau = [1 1 1];
SimuParam.tilt = [-45 0 45 5];                        
SimuParam.rfs = 0;

sim_param = SimuParam;

end