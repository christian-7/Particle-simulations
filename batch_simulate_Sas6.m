%% Generate library of rotated centrioles

clear; clc; close all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_of_structures   = 500;
current_dir         = ('Z:\Christian-Sieben\Centriole\data');
                                                                            % Set simulation parameters                                                           
labelling_eff   = 1;
nframes         = 20e3;
minNoise        = -200; 
maxNoise        = 200;
noise_mol       = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('Z:\Christian-Sieben\Centriole\data');
load('Sas6_extra_large.mat')

ptCloud = ptCloud_large;

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

fprintf('\n -------------------------------    Ready to start simulation ------------------------------- \n')

%% Plot a random example

close all
ID = randi([1 num_of_structures],1,1)
figure
scatter(sim_cent{ID, 1}(:,1),sim_cent{ID, 1}(:,2));
axis([-100 100 -100 100])

%% Run the simulator
tic
clc

fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,num_of_structures) '\n\n']);

parfor i = 1:length(sim_cent);
   
simCent_wNoise{i,1} = SMLM_simulator_batch(sim_cent{i,1}, nframes);

fprintf('\b|\n');

end

fprintf('\n -------------------------------    All Simulations Finished after %f  min ------------------------------- \n', toc/60)
toc
save('sim_Sas6_extra_large_noNoise_1.mat','simCent_wNoise');

%% Plot a random example along with GT

close all
ID = randi([1 num_of_structures],1,1)
figure
subplot(1,2,1)
scatter(sim_cent{ID, 1}(:,1),sim_cent{ID, 1}(:,2));
axis([-100 100 -100 100])

subplot(1,2,2)
scatter(simCent_wNoise{ID, 1}(:,1),simCent_wNoise{ID, 1}(:,2));
axis([-100 100 -100 100])

figure
scatter(simCent_wNoise{ID, 2}(:,1),simCent_wNoise{ID, 2}(:,2));
% axis([-100 100 -100 100])

%% Merging
tic
pos_list_new = {};

for j = 1:length(simCent_wNoise);

    pos_list      = [];
    pos_list(:,1) = simCent_wNoise{j,1}(:,1)+500; % get rid of negative values
    pos_list(:,2) = simCent_wNoise{j,1}(:,2)+500; % get rid of negative values
    pos_list(:,3) = simCent_wNoise{j,1}(:,3);
    pos_list(:,4) = simCent_wNoise{j,1}(:,4);
    pos_list = sortrows(pos_list,4);              % order by frames
    pos_list_new{j,1} = pos_list;             
    
end
    
simCent_wNoise{1,2} = merging(pos_list_new{1,1},15,50);

parfor i = 1:length(simCent_wNoise);
    

try
    
    simCent_wNoise{i,2} = merging(pos_list_new{i,1},15,50);
    
    catch
    
    simCent_wNoise{i,2} = [];

end
    
end


fprintf('\n -------------------------------    Data merged after %f  min ------------------------------- \n', toc/60)

%% Smart Merging
tic

simCent_wNoise{1,3} = smart_merging(pos_list_new{13,1},15,50);

parfor i = 1:length(simCent_wNoise);
    
try     
    simCent_wNoise{i,3} = smart_merging(pos_list_new{i,1},15,50);
catch
    simCent_wNoise{i,3} = [];
end

end

fprintf('\n -------------------------------    Data smart merged after %f  min ------------------------------- \n', toc/60)

%% Add Box around each Particle

data_column = 2; % locs =1 , merged =2;

% Define the bounding box 

xmin = -200;
ymin = -200;
xmax = 200;
ymax = 200;

interval = 10;

lowerLine = [];
lowerLine(:,1) = xmin:((max(xmax))/interval):max(xmax); 
lowerLine(1:end,2) = ymin; 

upperLine = [];
upperLine(:,1) = xmin:(max(xmax)/interval):max(xmax); 
upperLine(1:end,2) = max(ymax); 

leftLine = [];
leftLine(:,2) = ymin:((max(ymax))/interval):max(ymax); 
leftLine(1:end,1) = xmin; 

rightLine = [];
rightLine(:,2) = ymin:(max(ymax)/interval):max(ymax); 
rightLine(1:end,1) = max(xmax); 

addedBox = [lowerLine; upperLine; leftLine; rightLine];

figure
scatter(addedBox(:,1), addedBox(:,2));

% Add the box to each channel and separate the two channels

for i=1:length(simCent_wNoise);
    
     if isempty(simCent_wNoise{i,data_column});
     else
temp = [];         
temp(:,1) = simCent_wNoise{i,data_column}(:,1)-500;    
temp(:,2) = simCent_wNoise{i,data_column}(:,2)-500;
          
simCent_wNoise{i,data_column+1} = [temp; addedBox]; % Channel 1
     
     end
     
end


fprintf('-- Box added to each particle --');



%% Render each particle

data_column = 1; % locs =1 , merged =2;

pxlsize = 10;

Cent = simCent_wNoise;

current_dir = ('Z:\Christian-Sieben\Centriole\data');
[status, message, messageid] = rmdir('rendered', 's')
mkdir rendered
cd('rendered')

% Determine the box size form the largest particle

im_size = [];

for j=1:length(Cent);
    
    if isempty(Cent{j,data_column});
    else
        
idx = ~(isnan(Cent{j, data_column}(:,1))| isinf(Cent{j, data_column}(:,1))| isnan(Cent{j, data_column}(:,2))| isinf(Cent{j, data_column}(:,2)));

im_size(j,1) = round((max(Cent{j,data_column}(idx,1))-min(Cent{j,data_column}(idx,1)))/pxlsize);
im_size(j,2) = round((max(Cent{j,data_column}(idx,2))-min(Cent{j,data_column}(idx,2)))/pxlsize);
    
    end               
end

width  = max(max(im_size));
heigth = max(max(im_size));

for j = 1:length(Cent);
    
     if isempty(Cent{j,data_column});
     else
    
        idx = ~(isnan(Cent{j, data_column}(:,1))| isinf(Cent{j, data_column}(:,1))| isnan(Cent{j, data_column}(:,2))| isinf(Cent{j, data_column}(:,2)));

        rendered = hist3([Cent{j,data_column}(idx,2),Cent{j,data_column}(idx,1)],[max(max(im_size)) max(max(im_size))]);
        
        rendered_cropped = imcrop(rendered,[2 2 max(max(im_size))*0.9 max(max(im_size))*0.9]);
        
        empty  = zeros(round(max(max(im_size))*2),round(max(max(im_size))*2));
        
        size_DC     = size(empty);
        center_X    = round(size_DC(2)/2);
        center_Y    = round(size_DC(1)/2);

%         empty(round(center-heigth/2):(round(center-heigth/2))+heigth-1,round(center-width/2):(round(center-width/2))+width-1) = rendered_cropped;
        empty(round(center_Y-size(rendered_cropped,1)/2):(round(center_Y-size(rendered_cropped,1)/2))+size(rendered_cropped,1)-1,round(center_X-size(rendered_cropped,2)/2):(round(center_X-size(rendered_cropped,2)/2))+size(rendered_cropped,2)-1) = rendered_cropped;
        
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

%% Render Image without using the box
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_column         = 3;
image_column        = 6;
data_smart_merged   = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

current_dir = ('Z:\Christian-Sieben\Centriole\data');
[status, message, messageid] = rmdir('rendered', 's')
mkdir rendered
cd('rendered')

% Render each Particle

if data_smart_merged==1;
   image_column = data_column;
else
pxlsize = 10;

for i = 1:length(simCent_wNoise);
    
    if isempty(simCent_wNoise{i,data_column});
    else
    
    idx = ~(isnan(simCent_wNoise{i, data_column}(:,1))| isinf(simCent_wNoise{i, data_column}(:,1))| isnan(simCent_wNoise{i, data_column}(:,2))| isinf(simCent_wNoise{i, data_column}(:,2)));

heigth =  round((max(simCent_wNoise{i,data_column}(idx,1))-min(simCent_wNoise{i,data_column}(idx,1)))/pxlsize);
width  =  round((max(simCent_wNoise{i,data_column}(idx,2))-min(simCent_wNoise{i,data_column}(idx,2)))/pxlsize);

simCent_wNoise{i,image_column} = hist3([simCent_wNoise{i,data_column}(idx,2),simCent_wNoise{i,data_column}(idx,1)],[width heigth]); 
    
    end
end

end

% Insert each Particle into a black space

im_size = [];

for i = 1:length(simCent_wNoise);
    
im_size = vertcat(im_size,size(simCent_wNoise{i,image_column}));

end


for j = 1:length(simCent_wNoise);
    
       if isempty(simCent_wNoise{j,image_column});
       else

empty  = zeros(round(max(max(im_size))*2),round(max(max(im_size))*2));

heigth = size(simCent_wNoise{j,image_column},1);
width  = size(simCent_wNoise{j,image_column},2);
        
        size_DC     = size(empty);
        center      = round(size_DC(2)/2);
        

empty(round(center-heigth/2):(round(center-heigth/2))+heigth-1,round(center-width/2):(round(center-width/2))+width-1) = simCent_wNoise{j,image_column};

% empty(round(center_Y-size(rendered_cropped,1)/2):(round(center_Y-size(rendered_cropped,1)/2))+size(rendered_cropped,1)-1,round(center_X-size(rendered_cropped,2)/2):(round(center_X-size(rendered_cropped,2)/2))+size(rendered_cropped,2)-1) = rendered_cropped;
     

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
MIJ.run('Save', 'Tiff..., path=[Z:\\Christian-Sieben\\Centriole\\data\\new_woBox\\Montage_allParticles_Sas6_extra_large_1_noNoise_woBox_SmartMergedt.tif]');
MIJ.run('Close All')


% cd ..

fprintf('\n -------------------------------  Montage Saved  ------------------------------- \n')