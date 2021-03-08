%% Ideal Case

clear all; close all; clc

load('cam1_1.mat'); load('cam2_1.mat'); load('cam3_1.mat');
mass1_1 = trackMass(vidFrames1_1, 170, 430, 300, 400, 250);
mass2_1 = trackMass(vidFrames2_1, 100, 450, 200, 350, 250);
mass3_1 = trackMass(vidFrames3_1, 200, 350, 235, 500, 245);

% Put videos in sync
mass1_1 = trimmer(mass1_1);
mass2_1 = trimmer(mass2_1);
mass3_1 = trimmer(mass3_1);

case1 = [mass1_1'; mass2_1(1:length(mass1_1),:)'; mass3_1(1:length(mass1_1),:)'];

% Subtract the mean of each row from the data
[m,n] = size(case1); 
mean_of_row = mean(case1, 2); 
case1 = case1 - repmat(mean_of_row, 1, n); 

% Compute the SVD
[u,s,v] = svd(case1/sqrt(n-1));
DV = diag(s).^2; % diagonal variances
PC = u' * case1; % principal components

% plot
figure()
plot(1:6, DV/sum(DV), 'ro', 'Linewidth', 2);
title("Energy captured of each Diagonal Variance");
xlabel("Diagonal Variance"); 
ylabel("Energy");
set(gca, 'FontSize', 12)

figure()
subplot(2,1,1)
plot(1:length(mass1_1), case1(2,:),1:length(mass1_1), case1(1,:), 'Linewidth', 2)
title("Case 1: Displacement across Z axis and XY-plane");
xlabel("Time (frames)");
ylabel("Displacement (pixels)");
legend("Z", "XY")
set(gca, 'FontSize', 12)
subplot(2,1,2)
plot(1:length(mass1_1), PC(1,:),'Linewidth', 2)
title("Case 1: Displacement across Principal Component");
xlabel("Time (frames)"); 
ylabel("Displacement (pixels)");
legend("PC")
set(gca, 'FontSize', 12)

%% Noisy Case

clear all; close all; clc

load('cam1_2.mat'); load('cam2_2.mat'); load('cam3_2.mat');
mass1_2 = trackMass(vidFrames1_2, 170, 430, 300, 450, 250);
mass2_2 = trackMass(vidFrames2_2, 50, 475, 150, 450, 245);
mass3_2 = trackMass(vidFrames3_2, 200, 400, 210, 500, 245);

% Put videos in sync
mass1_2 = trimmer(mass1_2);
mass2_2 = trimmer(mass2_2);
mass3_2 = trimmer(mass3_2);

case2 = [mass1_2'; mass2_2(1:length(mass1_2),:)'; mass3_2(1:length(mass1_2),:)'];

% Subtract the mean of each row from the data
[m,n] = size(case2); 
mean_of_row = mean(case2, 2); 
case2 = case2 - repmat(mean_of_row, 1, n); 

% Compute the SVD
[u,s,v] = svd(case2/sqrt(n-1));
DV = diag(s).^2; % diagonal variances
PC = u' * case2; % principal components

% plot
figure()
plot(1:6, DV/sum(DV), 'ro', 'Linewidth', 2);
title("Energy captured of each Diagonal Variance");
xlabel("Diagonal Variance"); 
ylabel("Energy");
set(gca, 'FontSize', 12)

figure()
subplot(2,1,1)
plot(1:length(mass1_2), case2(2,:),1:length(mass1_2), case2(1,:), 'Linewidth', 2)
title("Case 2: Displacement across Z axis and XY-plane");
xlabel("Time (frames)");
ylabel("Displacement (pixels)");  
legend("Z", "XY")
set(gca, 'FontSize', 12)
subplot(2,1,2)
plot(1:length(mass1_2), PC(1,:), 1:length(mass1_2), PC(2,:), 'Linewidth', 2)
title("Case 2: Displacement across Principal Component");
xlabel("Time (frames)");
ylabel("Displacement (pixels)"); 
legend("PC1", "PC2")
set(gca, 'FontSize', 12)

%% Horizontal Displacement

clear all; close all; clc

load('cam1_3.mat'); load('cam2_3.mat'); load('cam3_3.mat');
mass1_3 = trackMass(vidFrames1_3, 225, 450, 275, 450, 250);
mass2_3 = trackMass(vidFrames2_3, 100, 450, 150, 425, 250);
mass3_3 = trackMass(vidFrames3_3, 150, 365, 210, 500, 245);

% Put videos in sync
mass1_3 = trimmer(mass1_3);
mass2_3 = trimmer(mass2_3);
mass3_3 = trimmer(mass3_3);

case3 = [mass1_3'; mass2_3(1:length(mass1_3),:)'; mass3_3(1:length(mass1_3),:)'];

% Subtract the mean of each row from the data
[m,n] = size(case3); 
mean_of_row = mean(case3, 2); 
case3 = case3 - repmat(mean_of_row, 1, n); 

% Compute the SVD
[u,s,v] = svd(case3/sqrt(n-1));
DV = diag(s).^2; % diagonal variances
PC = u' * case3; % principal components

% plot
figure()
plot(1:6, DV/sum(DV), 'ro', 'Linewidth', 2);
title("Energy captured of each Diagonal Variance");
xlabel("Diagonal Variance"); 
ylabel("Energy");
set(gca, 'FontSize', 12)

figure()
subplot(2,1,1)
plot(1:length(mass1_3), case3(2,:),1:length(mass1_3), case3(1,:), 'Linewidth', 2)
title("Case 3: Displacement across Z axis and XY-plane");
xlabel("Time (frames)");
ylabel("Displacement (pixels)");
legend("Z", "XY")
set(gca, 'FontSize', 12)
subplot(2,1,2)
plot(1:length(mass1_3), PC(1,:), 1:length(mass1_3), PC(2,:), 1:length(mass1_3), PC(3,:), 'Linewidth', 2)
title("Case 3: Displacement across Principal Component");
xlabel("Time (frames)");
ylabel("Displacement (pixels)");
legend("PC1", "PC2", "PC3")
set(gca, 'FontSize', 12)

%% Horizontal Displacement and Rotation

clear all; close all; clc

load('cam1_4.mat'); load('cam2_4.mat'); load('cam3_4.mat');
mass1_4 = trackMass(vidFrames1_4, 225, 450, 275, 450, 245);
mass2_4 = trackMass(vidFrames2_4, 100, 450, 150, 425, 245);
mass3_4 = trackMass(vidFrames3_4, 100, 300, 260, 500, 230);

% Put videos in sync
mass1_4 = trimmer(mass1_4);
mass2_4 = trimmer(mass2_4);
mass3_4 = trimmer(mass3_4);

case4 = [mass1_4(1:length(mass3_4),:)'; mass2_4(1:length(mass3_4),:)'; mass3_4'];

% Subtract the mean of each row from the data
[m,n] = size(case4); 
mean_of_row = mean(case4, 2); 
case4 = case4 - repmat(mean_of_row, 1, n); 

% Compute the SVD
[u,s,v] = svd(case4/sqrt(n-1));
DV = diag(s).^2; % diagonal variances
PC = u' * case4; % principal components

% plot
figure()
plot(1:6, DV/sum(DV), 'ro', 'Linewidth', 2);
title("Energy captured of each Diagonal Variance");
xlabel("Diagonal Variance"); 
ylabel("Energy");
set(gca, 'FontSize', 12)

figure()
subplot(2,1,1)
plot(1:length(mass3_4), case4(2,:),1:length(mass3_4), case4(1,:), 'Linewidth', 2)
title("Case 4: Displacement across Z axis and XY-plane");
xlabel("Time (frames)");
ylabel("Displacement (pixels)");
legend("Z", "XY")
set(gca, 'FontSize', 12)
subplot(2,1,2)
plot(1:length(mass3_4), PC(1,:), 1:length(mass3_4), PC(2,:), 1:length(mass3_4), PC(3,:), 'Linewidth', 2)
title("Case 4: Displacement across Principal Component");
xlabel("Time (frames)");
ylabel("Displacement (pixels)");
legend("PC1", "PC2", "PC3")
set(gca, 'FontSize', 12)


%% functions

function mass = trackMass(video, ymin, ymax, xmin, xmax, white)
    numOfFrames = size(video, 4); % number of frames in the given video
    mass = []; % x and y coordinates of the mass in each frame
    filter = zeros(480,640); 
    filter(ymin:1:ymax, xmin:1:xmax) = 1;
    for i = 1:numOfFrames
        I = rgb2gray(video(:,:,:,i)); 
        I = double(I).*filter; % zero out the values of pixels outside where we are interested in
        [Y, X] = ind2sub(size(I), find(I > white));
        mass = [mass; mean(X), mean(Y)]; % define the center of mass
    end
end

function sync_video = trimmer(raw_video)
    [M, I] = min(raw_video(1:30,2)); % start with the frame has the lowest y coordinate
    sync_video = raw_video(I:end,:);
end