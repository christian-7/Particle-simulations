%% Initialize two helices

clear; clc; close all

% generates z scale 

a  = 0: 20*pi/200 : 20*pi+5;
a2 = -5: 20*pi/200 : 20*pi;

% calculates x, y, and z

x = cos(a)*30;
y = sin(a)*30;
z = a/max(a)*250;

x2 = cos(a2)*30;
y2 = sin(a2)*30;
z2 = (a2+5)/max((a2+5))*250;

% Concatenate to molecule list

mol_list(:,1) = [x,x2]; 
mol_list(:,2) = [z2,z2]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

labelling_eff = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbr_of_labels = length(mol_list)*labelling_eff;

% Select random molecules from the list

r = randi([1,length(mol_list)],1,round(nbr_of_labels));

mol_list2 = [];
mol_list2 = mol_list(r,1:end);

% Plot the simulated points

figure('Position',[100 100 600 500])

width = 3;

subplot(2,2,1)
plot3(x,y,z,'LineWidth', width);hold on;
plot3(x2,y2,z2,'LineWidth', width)
grid on
title('Helix Structure - 3D view')
xlabel('[nm]')
ylabel('[nm]')
zlabel('[nm]')
axis([-80 80 -80 80 -30 280]);
box on

subplot(2,2,2)
plot3(x,y,z,'LineWidth', width);hold on;
plot3(x2,y2,z2,'LineWidth', width)
axis([-80 80 -80 80 -30 280]);
view(0,0)
grid on
title('Side View')
xlabel('[nm]')
ylabel('[nm]')
zlabel('[nm]')
box on

subplot(2,2,3)
scatter3(x,y,z,'*b');hold on;
scatter3(x2,y2,z2,'*r')
axis([-80 80 -80 80 -30 280]);
% view(-90,0)
grid on
title('Helix Points - 3D view')
xlabel('[nm]')
ylabel('[nm]')
zlabel('[nm]')
box on

subplot(2,2,4)
scatter(mol_list2(:,1),mol_list2(:,2));
grid on
axis([-80 80 -30 280]);
title(['NB position, LE = ', num2str(labelling_eff)]);
xlabel('[nm]')
ylabel('[nm]')
zlabel('[nm]')
box on
