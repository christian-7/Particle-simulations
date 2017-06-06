%% Generate library of rotated centrioles
clear; clc; close all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_of_structures   = 1000;
current_dir         = ('Z:\Christian-Sieben\Centriole\data');
                                                                            % Set simulation parameters                                                           
labelling_eff   = 0.0002;
nframes         = 20e3;
minNoise        = -400; 
maxNoise        = 400;
noise_mol       = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('Z:\Christian-Sieben\Centriole\data');
load('single_centriole.mat')

rng('shuffle');

a = -pi; b = pi;
sim_cent = {};

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NoiseX = (maxNoise-minNoise).*rand(1000,1) + minNoise; % short side, x position
NoiseY = (maxNoise-minNoise).*rand(1000,1) + minNoise; % short side, x position

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(sim_cent);
    
nbr_of_labels = length(sim_cent{i,1})*labelling_eff;
   
% Select random molecules from the list

r = randi([1,length(sim_cent{i,1})],1,round(nbr_of_labels));

% Add Noise molecules

NoiseMol(:,1) = NoiseX(randi([1 length(NoiseX)],noise_mol,1));
NoiseMol(:,2) = NoiseY(randi([1 length(NoiseY)],noise_mol,1));

sim_cent{i,2} = vertcat(sim_cent{i,1}(r,1:2),NoiseMol);             % this adds 10 noise molecules to mol_list

end

sim_cent(:,1) = [];    

fprintf('\n -------------------------------    Ready to start simulation ------------------------------- \n', toc/60)


%% Plot a random example

close all
ID = randi([1 num_of_structures],1,1)
figure
scatter(sim_cent{ID, 1}(:,1),sim_cent{ID, 1}(:,2));

%% Run the simulator
tic
clc

fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,num_of_structures) '\n\n']);

parfor i = 1:length(sim_cent);
   
simCent_wNoise{i,1} = SMLM_simulator_batch(sim_cent{i,1}, nframes);

fprintf('\b|\n');

% X = ['\n Finished ',num2str(i),' of ',num2str(length(sim_cent)),];
% clc;disp(X); 

end

fprintf('\n -------------------------------    All Simulations Finished after %f  min ------------------------------- \n', toc/60)
toc
save('simCent_noNoise_0002.mat','simCent_wNoise');
%% Render each particle

pxlsize = 10;
boxFactor =2;

Cent = simCent_wNoise;

current_dir = ('Z:\Christian-Sieben\Centriole\data');
[status, message, messageid] = rmdir('rendered', 's')
mkdir rendered
cd('rendered')

% Determine the box size form the largest particle

im_size = [];

for j=1:length(Cent);
    
    if isempty(Cent{j,1});
    else
        
idx = ~(isnan(Cent{j, 1}(:,1))| isinf(Cent{j, 1}(:,1))| isnan(Cent{j, 1}(:,2))| isinf(Cent{j, 1}(:,2)));

im_size(j,1) = round((max(Cent{j,1}(idx,1))-min(Cent{j,1}(idx,1)))/pxlsize);
im_size(j,2) = round((max(Cent{j,1}(idx,2))-min(Cent{j,1}(idx,2)))/pxlsize);
    
    end               
end

width  = max(max(im_size));
heigth = max(max(im_size));

for j = 1:length(Cent);
    
     if isempty(Cent{j,1});
     else
    
        idx = ~(isnan(Cent{j, 1}(:,1))| isinf(Cent{j, 1}(:,1))| isnan(Cent{j, 1}(:,2))| isinf(Cent{j, 1}(:,2)));

        rendered = hist3([Cent{j,1}(idx,2),Cent{j,1}(idx,1)],[max(max(im_size)) max(max(im_size))]);
        
        empty  = zeros(round(max(max(im_size))*boxFactor),round(max(max(im_size))*boxFactor));
        center = round(length(empty)/2);

        empty(round(center-heigth/2):(round(center-heigth/2))+heigth-1,round(center-width/2):(round(center-width/2))+width-1) = rendered;
        
%         emptyG = imgaussfilt(empty, 1);
   
name = ['image_10nm_per_pxl_' num2str(j) '.tiff'];

I32=[];
I32=uint32(empty);

t = Tiff(name,'w');

tagstruct.ImageLength     = size(I32,1);
tagstruct.ImageWidth      = size(I32,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
t.setTag(tagstruct)

t.write(I32);
t.close()
     end
end

fprintf('\n -------------------------------    Images Exported   ------------------------------- \n')

cd ..

%% Use ImageJ to create, blur and save the Montage

% addpath('D:\fiji-win64\Fiji.app\scripts');
% Miji;

save_path = ('Z:\Christian-Sieben\Centriole\data\rendered');

% cd('rendered');

string1 = ['Image Sequence...']; string2 = ['open=' num2str(save_path) ' file=image_10nm sort'];

MIJ.run(string1, string2);
MIJ.run('Gaussian Blur...','sigma=1 stack');

string1 = ['Make Montage...']; string2 = ['columns=' num2str(round(sqrt(length(simCent_wNoise)))) ' rows=' num2str(round(sqrt(length(simCent_wNoise))))  ' scale=1'];

MIJ.run(string1, string2);
MIJ.run('Enhance Contrast', 'saturated=0.35');
MIJ.run('Save', 'Tiff..., path=[Z:\\Christian-Sieben\\Centriole\\data\\Montage_allParticles_0002_noNoise_big_box.tif]');
MIJ.run('Close All')

cd ..

fprintf('\n -------------------------------  Montage Saved  ------------------------------- \n')