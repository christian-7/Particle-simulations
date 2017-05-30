%% Generate library of rotated centrioles

clear, clc

cd('Z:\Christian-Sieben\Centriole\data');
load('single_centriole.mat')

rng('shuffle');

a = -pi; b = pi;
sim_cent = {};
num_of_structures = 100;

rand_ang(:,1) = (b-a).*rand(num_of_structures,1) + a;
rand_ang(:,2) = (b-a).*rand(num_of_structures,1) + a;
rand_ang(:,3) = (b-a).*rand(num_of_structures,1) + a;

for i = 1:num_of_structures;

Rx = [1 0 0 0; ...
     0 cos(rand_ang(i,1)) -sin(rand_ang(i,1)) 0; ...
     0 sin(rand_ang(i,1)) cos(rand_ang(i,1)) 0; ...
     0 0 0 1];

tform = affine3d(Rx);
ptCloudOutx = pctransform(ptCloud,tform);
 
Ry = [cos(rand_ang(i,2)) 0 sin(rand_ang(i,2)) 0; ...
     0 1 0 0; ...
     -sin(rand_ang(i,2)) 0 cos(rand_ang(i,2)) 0; ...
     0 0 0 1];
 
tform = affine3d(Ry);
ptCloudOuty = pctransform(ptCloudOutx,tform);

Rz = [cos(rand_ang(i,3)) sin(rand_ang(i,3)) 0 0; ...
     -sin(rand_ang(i,3)) cos(rand_ang(i,3)) 0 0; ...
     0 0 1 0; ...
     0 0 0 1];

tform = affine3d(Rz);
ptCloudOut_final = pctransform(ptCloudOuty,tform);

% Substact the center of mass

sim_cent{i,1}(:,1) = ptCloudOut_final.Location(:,1)-sum(ptCloudOut_final.Location(:,1))/length(ptCloudOut_final.Location(:,1));
sim_cent{i,1}(:,2) = ptCloudOut_final.Location(:,2)-sum(ptCloudOut_final.Location(:,2))/length(ptCloudOut_final.Location(:,2));

clc; 
fprintf(['Done ' num2str(i) ' of ' num2str(num_of_structures)]);

end

cd('Z:\Christian-Sieben\Centriole\data');
save('sim_cent.mat','sim_cent');
clc; 
fprintf(['Library saved']);