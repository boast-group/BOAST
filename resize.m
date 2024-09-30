function res_vol3d = resize(vol3d, flag)

% Updated on 28/09/2024
% By Shokoufeh Golshani

% =========================================================================
% This function reduces the Field matrix size by a factor of 3 through
% averaging every three adjacent voxels
% =========================================================================

if flag == 0
    res_vol3d = vol3d;
else
    for r1 = 1:fix(size(vol3d, 3)/3)
        r1_vol3d(:,:,r1) = sum(vol3d(:,:,r1*3-2:r1*3), 3);
    end
 
    for r2 = 1:fix(size(vol3d,2)/3)
        r2_vol3d(:,r2,:) = sum(r1_vol3d(:,r2*3-2:r2*3,:), 2);
    end
 
    for r3 = 1:fix(size(vol3d,1)/3)
        res_vol3d(r3,:,:) = sum(r2_vol3d(r3*3-2:r3*3,:,:), 1)/27;
    end

end

end

