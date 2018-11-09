function r3_vol3d = myresize_3(vol3d,res)
% res = 0: no resize, res = 1: reducing size by factor 3

if res == 0
 r3_vol3d = vol3d;
else
 for r1 = 1:fix(size(vol3d,3)/3)
   r1_vol3d(:,:,r1) = sum(vol3d(:,:,r1*3-2:r1*3),3);
 end
 for r2 = 1:fix(size(vol3d,2)/3)
   r2_vol3d(:,r2,:) = sum(r1_vol3d(:,r2*3-2:r2*3,:),2);
 end
 for r3 = 1:fix(size(vol3d,1)/3)
   r3_vol3d(r3,:,:) = sum(r2_vol3d(r3*3-2:r3*3,:,:),1)/27;
 end
end

end
