%% Plot positions of Ab bound to the Sas6 cartwheel

% N-Terminus - towards the central hub
% C-Terminus - towards the MT barrel
% generate different radii to simulate different Ab positions

% Save the models as input for the simulation

clear, clc

r1 = 15; % nm
r2 = 45; % nm
% r2 = 25; % nm


  
center = [0,0];
MTx    = []; MTy = []; MTz = [];
t      = -pi:(2*pi/13):pi;
  
  for i = 1:10;
  
  MTx = vertcat(MTx,transpose(r1*cos(t)+center(1)));
  MTy = vertcat(MTy,transpose(r1*sin(t)+center(2)));
  MTz = vertcat(MTz,zeros(length(t),1)+i*15);
  
  end
  
  figure 
  scatter3(MTx,MTy,MTz,'blue','filled');hold on;
  
  ptCloud_small = pointCloud([MTx,MTy,MTz]);  MTx    = []; MTy = []; MTz = [];
  
  for i = 1:10;
  
  MTx = vertcat(MTx,transpose(r2*cos(t)+center(1)));
  MTy = vertcat(MTy,transpose(r2*sin(t)+center(2)));
  MTz = vertcat(MTz,zeros(length(t),1)+i*15);
  
  end
  
  scatter3(MTx,MTy,MTz,'red','filled');hold on;
  
  ptCloud_large = pointCloud([MTx,MTy,MTz]);
  
  figure
  pcshow(ptCloud_small);hold on;
  pcshow(ptCloud_large);
  
%% 

cd('Z:\Christian-Sieben\Centriole\data')    
save('Sas6_small.mat','ptCloud_small');
save('Sas6_extra_large.mat','ptCloud_large');