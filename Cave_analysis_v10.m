tic
clear; home
clf;
clc;
% Import data
window = 20; % Moving average window
threshold_percentile = 0.97; % Threshold to define the stalactites
threshold_flow_percentile = 0.93; % Threshold to define the possible flow area
aspect_ratio_fracture = 8;
aspect_ratio_matrix = 12;
eqv_area_solution_pipe = 100;
eqv_area_matrix = 75;
eqv_area_fracture = 100;
length_dia_ratio_soda = 4;
dia_soda = 1.5;

crop = LoadGrid('CWG(Mat)-Pillar-10mm-Sub Samp.SGEMS'); % Site 1 4mm
% crop2 = crop(1:4000, 1:1900);
%qz_crop = LoadGrid('qz_crop_4mm_site_2.SGEMS'); % Site 2 4mm
qz_crop = -crop; 

% qz_crop = LoadGrid('qz_crop_4mm_site_3.SGEMS'); % Site 3 4mm
% qz_crop = LoadGrid('qz_crop_4mm_site_1_leftportion.SGEMS'); % Site 1 leftportion 4mm
% qz_crop = LoadGrid('qz_crop_4mm_site_2_compare.SGEMS'); % Site 2 compare 4mm for validation
% qz_crop = qz_crop(:,51:end);

% qz_crop = LoadGrid('qz_crop_site3_slice1.SGEMS'); % Site 1
% qz_crop = LoadGrid('qz_crop_site_1_leftportion.SGEMS');
% qz_crop = LoadGrid('qz_crop_site_1_rightportion.SGEMS');
% qz_crop = LoadGrid('qz_crop_5mm_site_2_1.SGEMS'); % Site 2
% qz_crop = -qz_crop;
% qz_crop = LoadGrid('qz_crop_site_12_1.SGEMS'); % Site 3
% qz_crop = qz_crop(1:400,1:400);
% 
figure(1); clf;
imagesc(qz_crop);
colorbar
axis equal tight

%%
% % [data] = xlsread('/Golgotha_data_cropped_2.0.1.xlsx');
% % data = load( 'D:\Kashif\UNSW\PhD Research Work\Golgotha Cave\WA trip with LIDAR\Test data sets\Golgotha_site_3_slice_1.xyz');
% % data = load( 'D:\Kashif\UNSW\PhD Research Work\Golgotha Cave\WA trip with LIDAR\Test data sets\Golgotha_site_2_slice_5.xyz');
% % data = load( 'D:\Kashif\UNSW\PhD Research Work\Golgotha Cave\WA trip with LIDAR\Test data sets\Golgotha_site_12_3.xyz');
% 
% % data = load( 'D:\Kashif\UNSW\PhD Research Work\Golgotha Cave\WA trip with LIDAR\Test data sets\Golgotha_site_1_rightportion.xyz');
% data = load( 'D:\Kashif\UNSW\PhD Research Work\Golgotha Cave\WA trip with LIDAR\Test data sets\Golgotha_4mm_site_2_compare.xyz');
% % data = load( 'D:\Kashif\UNSW\PhD Research Work\Golgotha Cave\WA trip with LIDAR\Test data sets\Golgotha_small_test_slice_1.xyz');
% % data = load( 'D:\Kashif\UNSW\PhD Research Work\Golgotha Cave\WA trip with LIDAR\Test data sets\Golgotha_site_2_2.xyz');
% % data = load( 'D:\Kashif\UNSW\PhD Research Work\Golgotha Cave\WA trip with LIDAR\Test data sets\Golgotha_1mm_site_2.xyz');
% % data = load( 'D:\Kashif\UNSW\PhD Research Work\Golgotha Cave\WA trip with LIDAR\Test data sets\merged.txt');
% % 
% %% Read X, Y, Z values
% x = data(:,1); % X coordinates (meters)
% y = data(:,2); % Y coordinates (meters)
% z = data(:,3); % Z coordinates (meters)
% 
% %% Construct the interpolant
% F = TriScatteredInterp(x,y,z); % Performs similarly as griddata
% 
% min_x = min(x);
% max_x = max(x);
% int_x = 0.001; % Taking nodes of 1cm
% 
% min_y = min(y);
% max_y = max(y);
% int_y = 0.001; % Taking nodes of 1cm
% 
% %% Evaluate the interpolant at the locations (qx, qy). The corresponding value at these locations is qz
% ti_x = min_x:int_x:max_x; 
% ti_y = min_y:int_y:max_y; 
% [qx,qy] = meshgrid(ti_x,ti_y);
% qz = F(qx,qy);
% 
% % figure(1);clf;
% % mesh(qx,qy,qz);
% % hold on;
% % surf(qx,qy,qz);
% % % plot3(x,y,z,'.');
% % % scatter3(x,y,z);
% % axis tight
% % shading interp
% % colorbar
% % % view(60,-50)
% % view(-100,-10)
% % % view(0,-40)
% 
% %% View the whole image
% figure(1);clf;
% % subplot(1,3,1)
% imagesc(qz);
% colorbar
% axis equal tight
% 
% %% Rotate the image
% qz_1 = imrotate(qz,30); % Site 3
% % qz_1 = imrotate(qz,33); % Site 2
% % subplot(1,3,2)
% figure(2);clf;
% imagesc(qz_1);
% colorbar
% axis equal tight

%% Crop a rectangular section from the whole image
% qz_crop = qz_1(4501:11500,4501:12500); % for 1mm site 2; int_x, int_y = 0.001
% qz_crop = qz_1(1101:2900,1101:3100); % for 4mm site 2; int_x, int_y = 0.004
% qz_crop = qz(601:2800,401:1200); % for site 1 slice 3; int_x, int_y = 0.004
% qz_crop = qz(1201:5600,801:2400); % for site 1 slice 3; int_x, int_y = 0.002
% qz_crop = qz(71:870,61:750); % for site 2 slice 5; int_x, int_y = 0.01
% qz_crop = qz(176:2175,151:1875); % for site 2 slice 5; int_x, int_y = 0.004
% qz_crop = qz_1(401:850,251:750); % for site 12_2; int_x, int_y = 0.01
% qz_crop = qz_1(1001:2125,626:1875); % for site 12_2; int_x, int_y = 0.004

% qz_crop = qz(5810:27800,2060:20050); % for site 2_1; int_x, int_y = 0.01
% qz_crop = qz(1601:5300,401:4100); % for site 2_2; int_x, int_y = 0.001
% qz_crop = qz(1301:5300,401:4200); % for site 2_2; int_x, int_y = 0.005
% qz_crop = qz(71:870,61:750); % for site 2 slice 5; int_x, int_y = 0.01
% qz_crop = qz(201:700,301:1200); % for site 14_1; int_x, int_y = 0.01
% qz_crop = qz(71:390,21:450); % for site 15_1; int_x, int_y = 0.01
% qz_crop = qz(91:490,71:310); % for site 8_1; int_x, int_y = 0.01
% qz_crop = qz(81:380,101:330); % for site 12_1; int_x, int_y = 0.01
% qz_crop = qz_1(61:380,91:270); % for site_1_leftportion; int_x, int_y = 0.01
% qz_crop = qz_1(101:1100,251:700); % for site_1_leftportion; int_x, int_y = 0.004
% qz_crop = qz(311:760,123:522); % for site_1_rightportion; int_x, int_y = 0.01
% qz_crop = qz_1(61:390,71:250); % for site_1_compare; int_x, int_y = 0.004
% qz_crop = qz_1(201:440,181:520); % for site_1_compare_1; int_x, int_y = 0.004
% qz_crop = qz(2351:6850,751:4750); % for site_1_rightportion; int_x, int_y = 0.001
% qz_crop = qz_1(281:600,181:700); % for site_2_compare_1; int_x, int_y = 0.004
% qz_crop = qz_1(1121:2370,1021:2770); % for site_2_compare_1; int_x, int_y = 0.001

% zdir = [1 0 0];
% qz_1 = rotate(qz,zdir,25);
% qz_crop = qz_1(401:850,251:750); % for site 12_2; int_x, int_y = 0.01

% qz_crop = qz(12:116,38:88); % for slice 5; int_x, int_y = 0.01
% qz_crop = qz(121:1160,311:1110); % for slice 5; int_x, int_y = 0.001
% qz_crop = qz(61:300,61:200); % for slice 6; int_x, int_y = 0.01

% qz_crop = qz(51:300,11:210); % for site 1 slice 8; int_x, int_y = 0.01

% qz_crop = qz(20:229,19:218); % for site 2 slice 2; int_x, int_y = 0.01
% qz_crop = qz(201:2300,191:2185); % for site 2 slice 2; int_x, int_y = 0.001

% qz_crop = qz(200:1149,171:490); % for site 3 slice 1; int_x, int_y = 0.01
% qz_crop = qz(400:2298,342:980); % for site 3 slice 1; int_x, int_y = 0.005
% qz_crop = qz(2000:11490,1710:4900); % for site 3 slice 1; int_x, int_y = 0.001

% qz_crop = qz(51:220,21:200); % for small test slice 1; int_x, int_y = 0.001
% qz_crop = qz(6:25,4:21); % for small test slice 1; int_x, int_y = 0.01
% 
% figure(3);clf;
% % subplot(1,3,3)
% imagesc(qz_crop);
% colorbar
% axis equal tight
% % axis ij
% % xlabel('Distance(cm)')
% % ylabel('Distance(cm)')
% % WriteGrid(qz_crop,'qz_crop_site_1_rightportion.SGEMS', 'SGEMS');
% % WriteGrid(qz_crop,'qz_crop_1mm_site_2_compare.SGEMS', 'SGEMS');


%% Moving average to have a smooth surface
movavg=MovAvg(qz_crop,window);
% movavg = LoadGrid('Moving_average_site_2_1mm.SGEMS'); % Site 2
fprintf('\nMoving average window = %0.0f sq grid\n', window*2)
figure(3);clf
imagesc(movavg);
colorbar
axis equal tight
title('Moving average')
xlabel('Distance(cm)')
ylabel('Distance(cm)')
% WriteGrid(movavg,'Moving_average_4mm_site_1.SGEMS', 'SGEMS');
% WriteGrid(movavg,'Moving_average_4mm_site_2.SGEMS', 'SGEMS');
% WriteGrid(movavg,'Moving_average_4mm_site_3.SGEMS', 'SGEMS');

%% ceiling surface minus the moving average gives the topography
topo_initial = movavg-qz_crop;
% min(qz_crop(:))
% max(qz_crop(:))
% min(movavg(:))
% max(movavg(:))
% min(topo(:))
% max(topo(:))


%% Plot histogram of the anomalies - helps defining where to put the threshold
figure(4);clf;hold on
histogram(topo_initial(:), 100)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','k')
title('Histogram for Stalactites selection')
threshold = quantile(topo_initial(:),threshold_percentile);
vline(threshold,'r','97th percentile Threshold')
threshold_flow = quantile(topo_initial(:),threshold_flow_percentile);
vline(threshold_flow,'g','94th percentile Threshold')

% xmin = round(min(topo(:))*10)/10; xmax = round(max(topo(:))*10)/10;
% xlim([xmin xmax])
xlim([-0.2 0.3])
% [AX,H1] = plotyy(topo,'plot');
% set(AX,{'ycolor'},{'r'})
fprintf('Threshold Percentile = %0.2f Percent\n', threshold_percentile*100)
fprintf('Threshold = %0.4f\n', threshold)
fprintf('Threshold Percentile to identify Fracture flow = %0.2f Percent\n', threshold_flow_percentile*100)
fprintf('Threshold_flow = %0.4f\n', threshold_flow)

% %% anomalies in the surface (everything that departs from a smooth surface)
% figure(5);clf
% imagesc(topo_initial);
% colorbar
% colormap default
% axis equal tight
% title('Topography')
% xlabel('Distance(cm)')
% ylabel('Distance(cm)')
% 
% figure(6);clf
anomalies_initial=topo_initial>threshold; % Isolate all pixels with a deviation from the moving surface above a threshold
% imagesc(anomalies);
% axis equal tight
% colorbar
% colormap gray
% title('Anomalies/Location of Stalactites')
% xlabel('Distance(cm)')
% ylabel('Distance(cm)')


% %% Draw variogram of the ceiling topography anomalies to define the window size of moving average
% figure(55);clf;
% pts=Grid2Pts(anomalies,xx,yy,zeros(size(xx)));
% [meanh,gammah]=Variogram(pts,20,100,50,0);
% plot(meanh,gammah,'r','LineWidth',1)
% ylim([0 max(gammah)])
% axis equal square
% axis tight
% title('Variogram of the topography anomalies')
% xlabel('Distance')
% ylabel('Variance')


%% Defining the Locations of Stalactites
L_initial=bwlabel(anomalies_initial,8); % Identify the connected components. Each connected component in principle corresponds to one stal.

% imagesc(L_initial)

No_of_Stals_initial = max(L_initial(:));
% Density_of_stals = No_of_Stals/size(topo,1)*250/size(topo,2)*250;
% fprintf('Total Number of stalactites = %0.0f\n', No_of_Stals)
% fprintf('Ceiling dimension = %0.2f meter by %0.2f meter\n', size(anomalies,1)/250, size(anomalies,2)/250)
% fprintf('Density of stalactites = %0.0f per sq meter\n', Density_of_stals)

%% Use of regionprops function to find out 'Area','EquivDiameter','Orientation','Centroid'
STATS_1_initial = regionprops(logical(L_initial),'EquivDiameter');

EquivDiameter = (cat(1, STATS_1_initial.EquivDiameter)); 
Min_dia_exclude = min( cellfun(@(STATS_1_initial)min(STATS_1_initial(:)), {STATS_1_initial.EquivDiameter}) );
index = EquivDiameter > Min_dia_exclude;

exclude_loc = find((cat(1, STATS_1_initial.EquivDiameter)) == Min_dia_exclude);
I = zeros(size(exclude_loc));
J = zeros(size(exclude_loc));
anomalies = anomalies_initial;
for i=1:size(exclude_loc)
    [I(i),J(i)] = find(L_initial==exclude_loc(i));
    anomalies(I(i),J(i)) = 0;
end

L=bwlabel(anomalies,8);
STATS_1 = regionprops(logical(L),'Area','EquivDiameter','Orientation','Centroid','Image','BoundingBox');
imshow(L)
% Diameter = (cat(1, STATS_1.EquivDiameter))/2.5; % Stalactites diameters in cm units
% Area = (cat(1, STATS_1.Area))/2.5/2.5; % Stalactites diameters in cm units
Diameter = (cat(1, STATS_1.EquivDiameter))/10; % Stalactites diameters in cm units
Area = (cat(1, STATS_1.Area))/100; % Stalactites diameters in cm units

topo = topo_initial;
K = find(EquivDiameter == Min_dia_exclude);
for i=1:size(K,1)
    b = find(L_initial==K(i));
    topo(b) = 0;
%     B1(b1) = Aspect_Ratio_2(K1(i));
end

figure(5);clf;hold on
subplot(1,2,1)
imagesc(topo_initial);
colorbar
colormap default
axis equal tight
title('Initial Topography')
xlabel('Distance(cm)')
ylabel('Distance(cm)')
subplot(1,2,2)
imagesc(topo);
colorbar
colormap default
axis equal tight
title('Topography')
xlabel('Distance(cm)')
ylabel('Distance(cm)')

%%
%% Plot histogram of the anomalies - helps defining where to put the threshold
% figure(40);clf;hold on
% histogram(topo(:),100)
% % h = findobj(gca,'Type','patch');
% % set(h,'FaceColor','k')
% title('Histogram for Stalactites selection')
% threshold = quantile(topo(:),threshold_percentile);
% vline(threshold,'r','99th percentile Threshold')
% threshold_flow = quantile(topo(:),threshold_flow_percentile);
% vline(threshold_flow,'g','98th percentile Threshold')
% 
% % xmin = round(min(topo(:))*10)/10; xmax = round(max(topo(:))*10)/10;
% % xlim([xmin xmax])
% xlim([-0.2 0.3])
% % [AX,H1] = plotyy(topo,'plot');
% % set(AX,{'ycolor'},{'r'})
% fprintf('Threshold Percentile = %0.2f Percent\n', threshold_percentile*100)
% fprintf('Threshold = %0.4f\n', threshold)
% fprintf('Threshold Percentile to identify Fracture flow = %0.2f Percent\n', threshold_flow_percentile*100)


%%
figure(6);clf;hold on
% subplot(1,2,1)
% % anomalies_initial=topo_initial>threshold; % Isolate all pixels with a deviation from the moving surface above a threshold
% anomalies_initial_rot = imrotate(anomalies_initial,180);
% anomalies_initial_rot = flipdim(anomalies_initial_rot,1);
% imagesc(anomalies_initial_rot);
% axis equal tight
% colorbar
% colormap gray
% title('Initial Anomalies/Location of Stalactites')
% xlabel('Distance(cm)')
% ylabel('Distance(cm)')
% 
% subplot(1,2,2)
% anomalies=topo>threshold; % Isolate all pixels with a deviation from the moving surface above a threshold
anomalies_rot = imrotate(anomalies,180);
anomalies_rot = flipdim(anomalies_rot,1);
imagesc(anomalies_rot);
% imagesc(anomalies);
axis equal tight
colorbar
colormap gray
title('Anomalies/Location of Stalactites')
xlabel('Distance(cm)')
ylabel('Distance(cm)')

%%
% No_of_exact_Stals = size(Diameter,1);
No_of_Stals = max(L(:));

Density_of_stals = No_of_Stals/size(topo,1)*1000/size(topo,2)*1000;
fprintf('Total Number of stalactites = %0.0f\n', No_of_Stals)
fprintf('Ceiling dimension = %0.2f meter by %0.2f meter\n', size(anomalies,1)/1000, size(anomalies,2)/1000)
fprintf('Density of stalactites = %0.0f per sq meter\n', Density_of_stals)

% Area = Area(index);
Max_area = max(Area);
Min_area = min(Area);
% Min_area = (min( cellfun(@(STATS_1)max(STATS_1(:)), {STATS_1.Area}) ))/2.5/2.5;
fprintf('\nMaximum area of stalactite = %0.4f cm sq\n', Max_area)
fprintf('Minimum area of stalactite = %0.4f cm sq\n', Min_area)

Max_dia = max(Diameter);
Min_dia = min(Diameter);
% Max_dia = max( cellfun(@(STATS_1)max(STATS_1(:)), {STATS_1.EquivDiameter}) );
% Min_dia = min( cellfun(@(STATS_1)min(STATS_1(:)), {STATS_1.EquivDiameter}) );
fprintf('Maximum diameter of stalactite = %0.4f cm\n', Max_dia)
fprintf('Minimum diameter of stalactite = %0.4f cm\n', Min_dia)


%% Find out Stalactites lengths
% Equivlength = zeros(No_of_Stals);
% More appropriate way to find out the Stalactites lengths & Topographic elevations
ind = cell(1,No_of_Stals);
length_all = cell(1,No_of_Stals);
length = zeros(No_of_Stals,1);
elv_all = cell(1,No_of_Stals);
elv = zeros(No_of_Stals,1);
for i = 1:No_of_Stals
    ind{i} = find(L==i);
    length_all{i} = (topo(ind{i}));
    length(i) = max(length_all{i}); % Stalactites lengths in m units
    elv_all{i} = (qz_crop(ind{i}));
    elv(i) = max(elv_all{i})-length(i); % Topographic elevations in m units
end
length = abs((length - threshold)*100);
Sum_of_lengths = sum(length);
Avg_length = Sum_of_lengths/No_of_Stals; % Length in cm units
fprintf('\nSum of stalactites length = %0.2f meter\n', Sum_of_lengths)
fprintf('Average length of stalactites = %0.2f cm\n\n', Avg_length)
length_exclude_min = length;
elv_exclude_min = elv;


%% Find out the Aspect ratio, flow types
threshold_flow = quantile(topo_initial(:),threshold_flow_percentile);
fprintf('Threshold_flow = %0.4f\n', threshold_flow)
anomalies_flow=topo>threshold_flow;

% anomalies_flow=topo_initial>threshold_flow;
% anomalies_flow=topo>threshold;
L_flow=bwlabel(anomalies_flow,8);
% STATS_2_initial = regionprops(logical(L_flow),'EquivDiameter');
% EquivDiameter_2 = (cat(1, STATS_2_initial.EquivDiameter)); 
% % Min_dia_exclude_2 = min( cellfun(@(STATS_2)min(STATS_2(:)), {STATS_2.EquivDiameter}) );
% % index_2 = EquivDiameter_2 > Min_dia_exclude_2;
% % Area_2 = Area_2(index_2);
%  
% Min_dia_exclude_2 = min( cellfun(@(STATS_2)min(STATS_2(:)), {STATS_2_initial.EquivDiameter}) );
% index_2 = EquivDiameter_2 > Min_dia_exclude_2;
% 
% exclude_loc_2 = find((cat(1, STATS_2_initial.EquivDiameter)) == Min_dia_exclude_2);
% I = zeros(size(exclude_loc_2));
% J = zeros(size(exclude_loc_2));
% anomalies_2 = anomalies_initial;
% for i=1:size(exclude_loc)
%     [I(i),J(i)] = find(L_flow_initial==exclude_loc_2(i));
%     anomalies_2(I(i),J(i)) = 0;
% end
% 
% L_flow=bwlabel(anomalies_2,8);
STATS_2 = regionprops(logical(L_flow),'Area','EquivDiameter','Image','BoundingBox');
% Diameter_2 = (cat(1, STATS_2.EquivDiameter))/2.5; % Stalactites diameters in cm units
% Area_2 = (cat(1, STATS_2.Area))/2.5/2.5; % Stalactites diameters in cm units
Diameter_2 = (cat(1, STATS_2.EquivDiameter))/10; % Stalactites diameters in cm units
Area_2 = (cat(1, STATS_2.Area))/100; % Stalactites diameters in cm units

figure(7);clf;hold on
subplot(2,2,1)
imagesc(anomalies_flow)
axis equal tight
colormap gray
% colorbar
title('Locations of possible flow areas')
xlabel('Distance(cm)')
ylabel('Distance(cm)')

BoundingBox = cat(1, STATS_1.BoundingBox);
BoundingBox = BoundingBox(:,3:4);
AR_1 = zeros(max(L(:)),1);
for i = 1:max(L(:))
    AR_1(i) = BoundingBox(i,2) / BoundingBox(i,1);
end
Aspect_Ratio = AR_1/min(AR_1(:)); % Normalize the aspect ratio
% Aspect_Ratio = Aspect_Ratio(index);
Avg_aspect_Ratio = mean(Aspect_Ratio);


BoundingBox_2 = cat(1, STATS_2.BoundingBox);
BoundingBox_2 = BoundingBox_2(:,3:4);
AR_2 = zeros(max(L_flow(:)),1);
for i = 1:max(L_flow(:))
    AR_2(i) = BoundingBox_2(i,2) / BoundingBox_2(i,1);
end
Aspect_Ratio_2 = AR_2/min(AR_2(:)); % Normalize the aspect ratio
% Aspect_Ratio_2 = Aspect_Ratio_2(index_2);
Avg_aspect_Ratio_2 = mean(Aspect_Ratio_2);

% Define Locations of flow types
B1 = zeros(size(anomalies_flow));
% K1 = find(Aspect_Ratio_2 > aspect_ratio_fracture & Area_2 > eqv_area_fracture);
K1 = find(Aspect_Ratio_2 > aspect_ratio_fracture & Area_2 > eqv_area_fracture);
for i=1:size(K1,1)
    b1 = find(L_flow==K1(i));
    B1(b1) = 1;
%     B1(b1) = Aspect_Ratio_2(K1(i));
end
fprintf('Maximum aspect ratio of flow area = %0.1f \n', max(Aspect_Ratio_2))
fprintf('Minimum aspect ratio of flow area = %0.1f \n', min(Aspect_Ratio_2))
fprintf('Fracture flow locations are defined by both aspect ratio > %0.0f and Equivalent Area > %0.0f cm\n', ...
    aspect_ratio_fracture, eqv_area_fracture)
fprintf('Number of Fracture flow locations = %0.0f \n', size(K1,1))

subplot(2,2,2)
imagesc(B1)
axis equal tight
colormap gray
% colorbar
title('Locations of Pure Fracture flow')
xlabel('Distance(cm)')
ylabel('Distance(cm)')

B2 = zeros(size(anomalies));
K2 = find(Aspect_Ratio < aspect_ratio_matrix & Area < eqv_area_matrix);
for i=1:size(K2,1)
    b2 = find(L==K2(i));
%     B2(b2) = Aspect_Ratio(K2(i));
    B2(b2) = 1;
end
fprintf('Matrix flow locations are defined by both aspect ratio < %0.0f and Equivalent Area < %0.0f cm\n', ...
    aspect_ratio_matrix, eqv_area_matrix)
fprintf('Number of Matrix flow locations = %0.0f \n', size(K2,1))
fprintf('Maximum aspect ratio of stalactite = %0.1f \n', max(Aspect_Ratio))
fprintf('Minimum aspect ratio of stalactite = %0.1f \n', min(Aspect_Ratio))
fprintf('Average aspect ratio of stalactite = %0.1f \n', Avg_aspect_Ratio)

subplot(2,2,3)
B2_rot = imrotate(B2,180);
B2_rot = flipdim(B2_rot,1);
imagesc(B2_rot);
% imagesc(B2)
axis equal tight
colormap gray
% colorbar
title('Locations of Matrix flow')
xlabel('Distance(cm)')
ylabel('Distance(cm)')

B3 = zeros(size(anomalies_flow));
K3 = find(Aspect_Ratio_2 < aspect_ratio_fracture & Area_2 > eqv_area_solution_pipe);
for i=1:size(K3,1)
    b3 = find(L_flow==K3(i));
    B3(b3) = 1;
%     B3(b3) = Aspect_Ratio_2(K3(i));
end
fprintf('Combination of Pipe, Fracture and Martix flow locations are defined by both aspect ratio < %0.0f and Equivalent Area > %0.0f cm\n', ...
    aspect_ratio_fracture, eqv_area_solution_pipe)
fprintf('Number of Combined flow locations = %0.0f \n', size(K3,1))
subplot(2,2,4)
imagesc(B3)
axis equal tight
colormap gray
% colorbar
title('Locations of Combination of Pipe, Fracture and Martix flow')
xlabel('Distance(cm)')
ylabel('Distance(cm)')


%% Find out Stalactites lengths
% Equivlength = zeros(No_of_Stals);
% More appropriate way to find out the Stalactites lengths & Topographic elevations
No_of_flow = max(L_flow(:));
ind_flow = cell(1,No_of_flow);
length_all_flow = cell(1,No_of_flow);
length_flow = zeros(No_of_flow,1);
elv_all_flow = cell(1,No_of_flow);
elv_flow = zeros(No_of_flow,1);
for i = 1:No_of_flow
    ind{i} = find(L_flow==i);
    length_all_flow{i} = (topo(ind{i}));
    length_flow(i) = max(length_all_flow{i}); % Stalactites lengths in m units
    elv_all_flow{i} = (qz_crop(ind{i}));
    elv_flow(i) = max(elv_all_flow{i})-length_flow(i); % Topographic elevations in m units
end
length_flow = abs((length_flow - threshold_flow)*100);
% Sum_of_lengths = sum(length);
% Avg_length = Sum_of_lengths/No_of_Stals; % Length in cm units
% fprintf('\nSum of stalactites length = %0.2f meter\n', Sum_of_lengths)
% fprintf('Average length of stalactites = %0.2f cm\n\n', Avg_length)
% length_exclude_min = length;
% elv_exclude_min = elv;


%% Finding soda-straws
index_soda = Diameter_2 < dia_soda;
Diameter_soda = Diameter_2(index_soda); 
length_soda = length_flow(index_soda);
elv_soda = elv_flow(index_soda);
length_dia_ratio = length_flow./Diameter_2;

B4 = zeros(size(anomalies_flow));
% aa = length_dia_ratio > length_dia_ratio_soda;
% bb = Area_initial < dia_soda*pi/4;
% cc = Area_initial > min(Area_initial);

% K4 = find(length_dia_ratio > length_dia_ratio_soda & Area_2 < dia_soda*dia_soda*pi/4);
K4 = find(length_dia_ratio > length_dia_ratio_soda & Diameter_2 < dia_soda);
for i=1:size(K4,1)
    b4 = find(L_flow==K4(i));
%     B2(b2) = Aspect_Ratio(K2(i));
    B4(b4) = 1;
end
fprintf('Soda-straws are defined by both length diameter ratio > %0.0f and Diameter < %0.2f cm\n', ...
    length_dia_ratio_soda, dia_soda)
fprintf('Number of Soda-straw Stalactites location = %0.0f \n', size(K4,1))
figure(8);clf; hold on
B4_rot = imrotate(B4,180);
B4_rot = flipdim(B4_rot,1);
imagesc(B4_rot);
% imagesc(B4)
% [I,J] = find(B4 == 1);
% C = ones(size(I,1),1)*1;
% hold on
% scatter(J,I,1,C,'.','r')
axis equal tight
colormap gray
% colorbar
title('Locations of Soda-straws')
xlabel('Distance(cm)')
ylabel('Distance(cm)')

fprintf('Total Number of stalactites = %0.0f\n', No_of_Stals)
Tot_stals_check = sum(size(K1,1)+size(K2,1)+size(K3,1)+size(K4,1));
Density_of_stals = Tot_stals_check/size(topo,1)*1000/size(topo,2)*1000;
fprintf('Density of stalactites = %0.0f per sq meter\n', Density_of_stals)

Aspect_ratio_plot = Aspect_Ratio_2(K1);
length_plot = length_flow(K1);
Aspect_ratio_plot( size(K1,1)+1 : size(K1,1)+size(K3,1)) = Aspect_Ratio_2(K3);
length_plot( size(K1,1)+1 : size(K1,1)+size(K3,1)) = length_flow(K3);
Aspect_ratio_plot( size(K1,1)+size(K3,1)+1 : size(K1,1)+size(K3,1)+size(K4,1)) = Aspect_Ratio_2(K4);
length_plot( size(K1,1)+size(K3,1)+1 : size(K1,1)+size(K3,1)+size(K4,1)) = length_flow(K4);
% Aspect_ratio_plot( size(K1,1)+size(K4,1)+size(K3,1)+1 : size(K1,1)+size(K3,1)+size(K4,1)+size(K2,1)) = Aspect_Ratio(K2);
% length_plot( size(K1,1)+size(K4,1)+size(K3,1)+1 : size(K1,1)+size(K3,1)+size(K4,1)+size(K2,1)) = length(K2);

% % WriteGrid files for site 1
% WriteGrid(B1,'Locations of Fracture flow (Type 2) site 1_4mm_2.SGEMS', 'SGEMS'); %Locations of Pure Fracture flow
% WriteGrid(B2,'Locations of Matrix flow (Type 1) site 1_4mm_2.SGEMS', 'SGEMS'); %Locations of Matrix flow
% WriteGrid(B3,'Locations of Combined flow (Type 3) site 1_4mm_2.SGEMS', 'SGEMS'); %Locations of Combined flow
% WriteGrid(B4,'Locations of Soda straw stals site 1_4mm_2.SGEMS', 'SGEMS'); %Locations of Soda straw stals
% WriteGrid(length_plot,'length_site1.SGEMS', 'SGEMS');
% WriteGrid(Aspect_ratio_plot,'Aspect_Ratio_site1.SGEMS', 'SGEMS');

% % WriteGrid files for site 2
% WriteGrid(B1,'Locations of Fracture flow (Type 2) site 2_4mm_2.SGEMS', 'SGEMS'); %Locations of Pure Fracture flow
% WriteGrid(B2,'Locations of Matrix flow (Type 1) site 2_4mm_2.SGEMS', 'SGEMS'); %Locations of Matrix flow
% WriteGrid(B3,'Locations of Combined flow (Type 3) site 2_4mm_2.SGEMS', 'SGEMS'); %Locations of Combined flow
% WriteGrid(B4,'Locations of Soda straw stals site 2_4mm_2.SGEMS', 'SGEMS'); %Locations of Soda straw stals
% WriteGrid(length_plot,'length_site2.SGEMS', 'SGEMS');
% WriteGrid(Aspect_ratio_plot,'Aspect_Ratio_site2.SGEMS', 'SGEMS');

% % WriteGrid files for site 3
% WriteGrid(B1,'Locations of Fracture flow (Type 2) site 3_4mm_2.SGEMS', 'SGEMS'); %Locations of Pure Fracture flow
% WriteGrid(B2,'Locations of Matrix flow (Type 1) site 3_4mm_2.SGEMS', 'SGEMS'); %Locations of Matrix flow
% WriteGrid(B3,'Locations of Combined flow (Type 3) site 3_4mm_2.SGEMS', 'SGEMS'); %Locations of Combined flow
% WriteGrid(B4,'Locations of Soda straw stals site 3_4mm_2.SGEMS', 'SGEMS'); %Locations of Soda straw stals
% WriteGrid(length_plot,'length_site3.SGEMS', 'SGEMS');
% WriteGrid(Aspect_ratio_plot,'Aspect_Ratio_site3.SGEMS', 'SGEMS');


%% Plotting Stalactites length vs diameter
% remove the stals having 1 grid of sizes which makes the log-log plot uneven 
% this is only for log-log plot

figure(9);clf;hold on
param_1 = zeros(size(length_exclude_min,1),2);
param_1(:,1) = log(length_exclude_min);
param_1(:,2) = log(Diameter);
% ksr(log(length),log(Diameter));
% ksr(length,Diameter);
% [f,xi] = ksdensity(values);
% plot(xi,f);

[pdf_1,X_1,Y_1] = Ksmooth(param_1,(std(param_1(:,1)))/3,(std(param_1(:,2)))/3,[min(param_1(:,1)) 0.01 max(param_1(:,1))],[min(param_1(:,2)) 0.01 max(param_1(:,2))]);
pdf_1 = pdf_1/sum(pdf_1(:));
surf(X_1,Y_1,pdf_1)
% mesh(X,Y,pdf,'EdgeColor','black')
shading flat
% axis equal square
axis equal tight
colormap default
colorbar
title('Plot for Stalactites with Soda-straws')
% xlim([-1.2 3.8])
% ylim([0.5 3])
xlabel('log length [log(cm)]')
ylabel('log diameter [log(cm)]')

% correlation between stalactites diameter and length
r = (corr(log(length_exclude_min),log(Diameter)));
fprintf('Correlation (r) between stals diameter and length = %0.2f\n', r)
% rr2 = (corr(length,Diameter))^2;
% slope = param_1(:,1)./param_1(:,2);
% avg_slope = mean(slope)

slope = gradient(param_1);
slope = abs(slope);
avg_slope = median(slope(:,1));
interquartile_range = iqr(slope(:,1));

% avg_slope = mean(slope(:,1))
% Avg_dia = sum(Diameter)/No_of_Stals;
% Percentile_5 = prctile(slope(:,1),5)
% upper_limit = avg_slope + Percentile_5;
% lower_limit = avg_slope - Percentile_5;

% mx=mean(slope(:,1));
% sx=std(slope(:,1));
% upper_limit = mx+1.96*sx;
% lower_limit = mx-1.96*sx;
% if lower_limit <= 0
%     lower_limit = 0;
% end
% figure(20);clf;
% plot(1:size(slope,1),slope(:,1),'o',1:size(slope,1),upper_limit,1:size(slope,1),lower_limit)


% Percentile_95 = prctile(slope(:,1),95)
% Percentile_95_1 = quantile(slope(:,1),0.95)

index_soda = Diameter > dia_soda;
% Diameter_soda = Diameter(index_soda); 
% length_soda = length(index_soda);
% elv_soda = elv_flow(index_soda);
Diameter_exclude_soda = Diameter(index_soda);
length_exclude_soda = length_exclude_min(index_soda);
elv_exclude_soda = elv_exclude_min(index_soda);

figure(10);clf;hold on
param_1 = zeros(size(length_exclude_soda,1),2);
param_1(:,1) = log(length_exclude_soda);
param_1(:,2) = log(Diameter_exclude_soda);

[pdf_1,X_1,Y_1] = Ksmooth(param_1,(std(param_1(:,1)))/3,(std(param_1(:,2)))/3,[min(param_1(:,1)) 0.01 max(param_1(:,1))],[min(param_1(:,2)) 0.01 max(param_1(:,2))]);
pdf_1 = pdf_1/sum(pdf_1(:));
surf(X_1,Y_1,pdf_1)
shading flat
% axis equal square
axis tight
colormap default
colorbar
title('Plot for Stalactites without Soda-straws')
% xlim([-1.2 4])
% ylim([0.1 2.8])
xlabel('log length [log(cm)]')
ylabel('log diameter [log(cm)]')

% correlation between stalactites diameter and length
r = (corr(log(length_exclude_soda),log(Diameter_exclude_soda)));
fprintf('Correlation (r) between stals diameter and length = %0.2f\n', r)
% rr2 = (corr(length,Diameter))^2;

slope_exclude_soda = gradient(param_1);
slope_exclude_soda = abs(slope_exclude_soda);
avg_slope_exclude_soda = median(slope_exclude_soda(:,1));
interquartile_range_exclude_soda = iqr(slope_exclude_soda(:,1));


%% Plotting Stalactites length vs topography elevation
figure(11);clf;hold on
param_2 = zeros(size(length_exclude_min,1),2);
param_2(:,1) = length_exclude_min;
param_2(:,2) = elv_exclude_min;
[pdf_2,X_2,Y_2] = Ksmooth(param_2,(std(param_2(:,1)))/3,(std(param_2(:,2)))/3,[min(param_2(:,1)) 0.01 max(param_2(:,1))],[min(param_2(:,2)) 0.01 max(param_2(:,2))]);
pdf_2 = pdf_2/sum(pdf_2(:));

surf(X_2,Y_2,pdf_2)
shading flat
% plot(length,elv,'.')

% plot([min(length),max(length)],[min(elv),max(elv)],'g')
% param_3(:,1) = log(length);
% param_4(:,2) = log(elv);
% 
% [pdf_1,X_1,Y_1] = Ksmooth(length,(std(length))/7,(std(elv))/7,[min(length(:,1)) 0.01 max(length)],[min(elv) 0.01 max(elv)]);
% scatter(X(:),Y(:),5,pdf(:))
% surf(X_1,Y_1,pdf_1)
% mesh(X,Y,pdf,'EdgeColor','black')

colormap default
% axis equal square
axis tight
colorbar
% title('Stals Length vs Elevation')
% xlim([0 30])
xlabel('Length (cm)')
ylabel('Elevation (m)')



%% Stalactite length vs. Aspect Ratio

figure(12); clf; hold on
scatter(length, Aspect_Ratio, 'red', 'filled')
xlim([0 100])
ylim([0 40])
title('site 1')
xlabel('Stalactite length (cm)')
ylabel('Aspect Ratio')
axis tight

%% Combination of flow types on single graph

ind_soda_s1 = find(B4==1);
B1(ind_soda_s1) = 4;
Percent_soda_straw_Site1 = size(ind_soda_s1,1)/size(B1,1)/size(B1,2)*100;
ind = find(B2==1);
B1(ind) = 2;
Percent_Flow2_Site1 = size(ind,1)/size(B1,1)/size(B1,2)*100; %matrix
ind = find(B3==1); 
B1(ind) = 3;
Percent_Flow3_Site1 = size(ind,1)/size(B1,1)/size(B1,2)*100; %fracture
ind = find(B1==1);
Percent_Flow1_Site1 = size(ind,1)/size(B1,1)/size(B1,2)*100; %combination
Total_flow_area_site1 = Percent_soda_straw_Site1 + Percent_Flow1_Site1 + Percent_Flow2_Site1 + Percent_Flow3_Site1
Proportion_type1_site1 = Percent_Flow1_Site1 / Total_flow_area_site1 * 100
Proportion_type2_site1 = Percent_Flow2_Site1 / Total_flow_area_site1 * 100
Proportion_type3_site1 = Percent_Flow3_Site1 / Total_flow_area_site1 * 100
Proportion_type4_site1 = Percent_soda_straw_Site1 / Total_flow_area_site1 * 100
No_flow_area_of_ceiling_site1 = (100 - Total_flow_area_site1) 



ylabel('Distance(cm)')
figure(13); clf; 
% subplot(1,3,1)
B1_rot = imrotate(B1,-90);
B1_rot = flipdim(B1_rot,1);
imagesc(B1_rot);
[I,J] = find(B1 == 4);
C = ones(size(I,1),1)*4;
hold on
scatter(size(B1,1)-I,size(B1,2)-J,.005,C,'*','r')
axis equal tight
cmap = jet(5);
cmap(1,:)=[0,0,0];
colormap(cmap)

hold on
L = line(ones(5),ones(5), 'LineWidth',2);               % generate line
set(L,{'color'},mat2cell(cmap,ones(1,5),3));            % set the colors according to cmap
%legend('Soda straw stals','No flow','Flow type 1','Flow type 2','Flow type 3','Location','SouthWest')

title('Locations of flow types for site 1')
xlabel('Distance(cm)')

% % figure(13);clf
% subplot(1,3,2)
% imagesc(B1_2)
% [I,J] = find(B1_2 == 4);
% C = ones(size(I,1),1)*4;
% hold on
% scatter(J,I,3,C,'*','r')
% axis equal tight
% colormap(cmap)
% title('Locations of flow types for site 2')
% xlabel('Distance(cm)')
% ylabel('Distance(cm)')
% 
% % figure(14);clf
% subplot(1,3,3)
% imagesc(B1_3)
% [I,J] = find(B1_3 == 4);
% C = ones(size(I,1),1)*4;
% hold on
% scatter(J,I,2,C,'*','r')
% axis equal tight
% colormap(cmap)
% title('Locations of flow types for site 3')
% xlabel('Distance(cm)')
% ylabel('Distance(cm)')

ylabel('Distance(cm)')
figure(14); clf; 
% subplot(1,3,1)
B1_rot = imrotate(B1,-90);
B1_rot = flipdim(B1_rot,1);
imagesc(B1_rot);
[I,J] = find(B1 == 4);
C = ones(size(I,1),1)*4;
hold on
scatter(size(B1,1)-I,size(B1,2)-J,.006,C,'*','r')
axis equal tight
cmap = jet(5);
cmap(1,:)=[0,0,0];
colormap(cmap)

hold on
L = line(ones(5),ones(5), 'LineWidth',2);               % generate line
set(L,{'color'},mat2cell(cmap,ones(1,5),3));            % set the colors according to cmap
%legend('Soda straw stals','No flow','Flow type 1','Flow type 2','Flow type 3','Location','SouthWest')

title('Locations of flow types for site 1')
xlabel('Distance(cm)')

ylabel('Distance(cm)')
figure(15); clf; 
% subplot(1,3,1)
B1_rot = imrotate(B1,-90);
B1_rot = flipdim(B1_rot,1);
imagesc(B1_rot);
[I,J] = find(B1 == 4);
C = ones(size(I,1),1)*4;
hold on
scatter(size(B1,1)-I,size(B1,2)-J,.007,C,'*','r')
axis equal tight
cmap = jet(5);
cmap(1,:)=[0,0,0];
colormap(cmap)

hold on
L = line(ones(5),ones(5), 'LineWidth',2);               % generate line
set(L,{'color'},mat2cell(cmap,ones(1,5),3));            % set the colors according to cmap
%legend('Soda straw stals','No flow','Flow type 1','Flow type 2','Flow type 3','Location','SouthWest')

title('Locations of flow types for site 1')
xlabel('Distance(cm)')

toc
beep