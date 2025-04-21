tic
clear;
clc;

% open excel file
[data]=readmatrix('PC-3-Sub Samp-10.txt');

% read X, Y, Z values
x = data(:,1); % X coordinates (meters)
y = data(:,2); % Y coordinates (meters)
z = data(:,3); % Z coordinates (meters)

 % Construct the interpolant
F = scatteredInterpolant(x,y,z,'natural','none'); % Performs similarly as griddata
min_x = min(x);
max_x = max(x);
int_x = 0.001; % Taking nodes of 1mm

min_y = min(y);
max_y = max(y);
int_y = 0.001; % Taking nodes of 1mm

%% Evaluate the interpolant at the locations (qx, qy). The corresponding value at these locations is qz
ti_x = min_x:int_x:max_x; 
ti_y = min_y:int_y:max_y; 
[qx,qy] = meshgrid(ti_x,ti_y);
qz = F(qx,qy); % 

figure(20);clf;
mesh(qx,qy,qz);
hold on;
surf(qx,qy,qz);
%plot3(x,y,z,'.');
% scatter3(x,y,z);
axis tight
shading interp
colorbar
colormap("hot")
%% View the whole image
figure(21);clf;
% subplot(1,3,1)
imagesc(qz);
colorbar
colormap("hot")

axis equal tight

%WriteGrid(qz, 'PC-1-Near-None.SGEMS', 'SGEMS');
toc
beep