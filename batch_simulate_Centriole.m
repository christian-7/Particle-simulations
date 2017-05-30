clear; clc; close all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('Z:\Christian-Sieben\Centriole\data');                           % Load single Centriole GT
load('sim_cent.mat')

save_path = 'Z:\Christian-Sieben\Centriole\data\rendered';          % Path for the images

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                    % Set simulation parameters
                                                               
labelling_eff   = 0.001;
nframes         = 20e3;
minNoise        = -400; 
maxNoise        = 400;
noise_mol       = 15;

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

fprintf(' -------------------------------    Ready to start simulation -------------------------------', toc/60)

%% Run the simulator
tic
parfor i = 1:length(sim_cent);
   
simCent_wNoise{i,1} = SMLM_simulator_batch(sim_cent{i,2}, nframes);

X = [' Finished ',num2str(i),' of ',num2str(length(sim_cent)),];
clc;disp(X); 

end

fprintf(' -------------------------------    All Simulations Finished after %f  min -------------------------------', toc/60)
toc

%% Render each particle

pxlsize = 10;

Cent = simCent_wNoise;

cd('Z:\Christian-Sieben\Centriole\data\rendered');

% Determine the box size form the largest particle

im_size = [];

for j=1:length(Cent);


im_size(i,1) = round((max(Cent{j,1}(:,2))-min(Cent{j,1}(:,2)))/pxlsize);
im_size(i,2) = round((max(Cent{j,1}(:,1))-min(Cent{j,1}(:,1)))/pxlsize);
        
          
end

width  = max(max(im_size));
heigth = max(max(im_size));

for j = 1:length(Cent);
    
        idx = ~(isnan(Cent{j, 1}(:,1))| isinf(Cent{j, 1}(:,1))| isnan(Cent{j, 1}(:,2))| isinf(Cent{j, 1}(:,2)));

        rendered = hist3([Cent{j,1}(idx,2),Cent{j,1}(idx,1)],[max(max(im_size)) max(max(im_size))]);
        
        empty  = zeros(round(max(max(im_size))*1.5),round(max(max(im_size))*1.5));
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

disp(' -------------------------------    Images Exported   -------------------------------')

%% Use ImageJ to create, blur and save the Montage

% addpath('D:\fiji-win64\Fiji.app\scripts');
% Miji;

cd(save_path);

string1 = ['Image Sequence...']; string2 = ['open=' num2str(save_path) ' file=image_10nm sort'];

MIJ.run(string1, string2);
MIJ.run('Gaussian Blur...','sigma=1 stack');

string1 = ['Make Montage...']; string2 = ['columns=' num2str(round(sqrt(length(simCent_wNoise)))) ' rows=' num2str(round(sqrt(length(simCent_wNoise))))  ' scale=1'];

MIJ.run(string1, string2);
MIJ.run('Enhance Contrast', 'saturated=0.35');
MIJ.run('Save', 'Tiff..., path=[Z:\\Christian-Sieben\\Centriole\\data\\Montage_allParticles_cropped.tif]');
MIJ.run('Close All')

cd ..

disp(' -------------------------------  Montage Saved  -------------------------------')