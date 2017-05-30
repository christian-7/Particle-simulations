%% Generate single centriole (surface)

clear;clc;close all

list    = zeros(30,2);
diam    = 20;                   % Diameter of one tubulin cylinder
r       = 90;                   % Diameter of the centriole
inc1    = 0.98;
inc2    = -0.9;
center  = [0,0];


t=-pi:(2*pi/9):pi;
x=r*cos(t)+center(1);
y=r*sin(t)+center(2);
scatter(x,y);hold on;

r=85+diam;
x2=r*cos(t+inc1)+center(1);
y2=r*sin(t+inc1)+center(2);
scatter(x2,y2);

r=110+diam;
x3=r*cos(t+inc2)+center(1);
y3=r*sin(t+inc2)+center(2);
scatter(x3,y3);

list(1:10,1)=x;
list(11:20,1)=x2;
list(21:30,1)=x3;
list(1:10,2)=y;
list(11:20,2)=y2;
list(21:30,2)=y3;


for i=1:30;
[X,Y,Z]=cylinder(diam);
surf1 = surf(X+list(i,1),Y+list(i,2),Z*450); hold on;
shading interp

% 
% name = ['object' num2str(i) '.obj'];
% saveobjmesh(name,X,Y,Z);

end
  
%% Generate a single MT 
  
  clear,clc
  
  r      = 13;
  center = [0,0];
  MTx    = []; MTy = []; MTz = [];
  t      = -pi:(2*pi/13):pi;
  
  for i = 1:100;
  
  MTx = vertcat(MTx,transpose(r*cos(t)+center(1)));
  MTy = vertcat(MTy,transpose(r*sin(t)+center(2)));
  MTz = vertcat(MTz,zeros(length(t),1)+i*10);
  
  end
 
  scatter3(MTx,MTy,MTz)
  
%   forMarcel = [];
%   forMarcel(:,1) = MTy;
%   forMarcel(:,2) = MTz;
%   
%   cd('/Users/christian/Documents/Arbeit/MatLab/centriole_sim')
%   dlmwrite('singleMT.csv',forMarcel);
  
%% Generate single MT centriole barrel
 
clear,clc

% Generate list of MT centers

list    = zeros(30,2);
diam    = 20;                   % Diameter of one tubulin cylinder
r       = 90;                   % Diameter of the centriole
inc1    = 0.98;
inc2    = -0.9;
center  = [0,0];

t=-pi:(2*pi/9):pi;
x=r*cos(t)+center(1);
y=r*sin(t)+center(2);

r=85+diam;
x2=r*cos(t+inc1)+center(1);
y2=r*sin(t+inc1)+center(2);

r=110+diam;
x3=r*cos(t+inc2)+center(1);
y3=r*sin(t+inc2)+center(2);

list(1:10,1)=x;
list(11:20,1)=x2;
list(21:30,1)=x3;
list(1:10,2)=y;
list(11:20,2)=y2;
list(21:30,2)=y3;

scatter(list(:,1),list(:,2),'*');

% Generate tube around each MT center

allMTx = []; allMTy = []; allMTz = [];MTx = []; MTy = []; MTz = [];
r2 = 13;
t  = -pi:(2*pi/13):pi;

for i = 1:length(list);
  
  for j = 1:50;       
  
  MTx = vertcat(MTx,transpose(r2*cos(t)+list(i,1)));
  MTy = vertcat(MTy,transpose(r2*sin(t)+list(i,2)));
  MTz = vertcat(MTz,zeros(length(t),1)+j*10);
  
  end

  allMTx = vertcat(allMTx,MTx);
  allMTy = vertcat(allMTy,MTy);
  allMTz = vertcat(allMTz,MTz);
  
end
  
scatter3(allMTx,allMTy,allMTz,'*');

ptCloud = pointCloud([allMTx,allMTy,allMTz]);
pcshow(ptCloud);
 
%% Generate library of rotated centrioles

rng('shuffle');

a = -pi; b = pi;
sim_cent = {};
num_of_structures = 1;

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

sim_cent{i,1} = ptCloudOut_final.Location(:,1:3);

end
%% Plot single rotated example

subplot(1,2,1)
pcshow(ptCloud);
subplot(1,2,2)
pcshow(ptCloudOut_final);

figure
scatter3(ptCloudOut_final.Location(:,1),ptCloudOut_final.Location(:,2),ptCloudOut_final.Location(:,3),'.');
