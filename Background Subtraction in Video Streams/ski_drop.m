clear all; close all; clc

ski_drop = VideoReader('ski_drop.mp4');
ski_drop_frames = read(ski_drop);

%%

dt = 1 / ski_drop.FrameRate; % number of seconds per video frame
t = 0:dt:ski_drop.Duration;
nFrames = ski_drop.NumFrames;

for i = 1:nFrames
    I = rgb2gray(ski_drop_frames(:,:,:,i));
    I = im2double(I);
    X(:,i) = reshape(I,[],1);;
end

%% Create DMD matrices

X1 = X(:,1:end-1);
X2 = X(:,2:end);

%% SVD of X1

[U, S, V] = svd(X1,'econ');

%% Plot the singular value spectrum

plot(diag(S),'ko','Linewidth',2)
title('Singular Value Spectrum')
set(gca,'Fontsize', 14, 'Xlim', [0 378])
xlabel('singular indices')
ylabel('singular values')

%% Apply DMD

r = 20;
U_r = U(:, 1:r); % truncate to rank-r
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);

Stilde = U_r' * X2 * V_r / S_r;
[eV, D] = eig(Stilde); % compute eigenvalues + eigenvectors
mu = diag(D); % extract eigenvalues
omega = log(mu)/dt;

plot(real(omega), imag(omega), 'ro', 'MarkerFaceColor', 'r')
title('Omega')
set(gca, 'Fontsize', 14 , 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
xlabel('Re')
ylabel('Im')

%% 

Phi = U_r * eV;
b = Phi \ X1(:,1);
ind = find(abs(omega) < 0.001);

%% reconstruction

U_modes = [];
for iter = 1:length(t)
   U_modes(:, iter) = b .* exp(omega * t(iter)); 
end
U_dmd = Phi(:,14) * U_modes(14,:);

%%
video_background = U_dmd;
video_foreground = X - abs(U_dmd);
video_foreground = (video_foreground - min(video_foreground)) ./ (max(video_foreground) - min(video_foreground));

%% sample original video

sample_video_original = reshape(X(:, 200), 540, 960); 
imshow(sample_video_original);

%% sample background

sample_video_background = reshape(video_background(:,200), 540, 960);
imshow((sample_video_background));
%% sample foreground

sample_video_foreground = reshape(video_foreground(:,200), 540, 960);
imshow(sample_video_foreground);