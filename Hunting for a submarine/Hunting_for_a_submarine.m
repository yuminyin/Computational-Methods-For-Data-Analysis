clear all; close all; clc
load subdata.mat

L = 10; % spatial domain
n = 64; % Fourier modes

x2 = linspace(-L,L,n+1); x = x2(1:n); y = x; z = x;
k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; ks = fftshift(k);

[X,Y,Z]=meshgrid(x,y,z); % returns 3-D grid coordinates defined by the vectors x, y, and z
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

% Average of the spectrum
ave = zeros(n,n,n);
for i = 1:49
    ave = ave + fftn(reshape(subdata(:,i),n,n,n)); % reshape and apply the Fourier transform
end
ave = abs(fftshift(ave))/49; % shift the Fourier transform data and take an average over realizations

% Determine the center frequency and its location in the frequency domain
maxAve = max(abs(ave),[],'all');
[x_index,y_index,z_index] = ind2sub(size(ave), find(abs(ave) == maxAve));
maxFrequencyAt = [Kx(x_index,y_index,z_index),Ky(x_index,y_index,z_index),Kz(x_index,y_index,z_index)];
x0 = Kx(x_index,y_index,z_index)
y0 = Ky(x_index,y_index,z_index)
z0 = Kz(x_index,y_index,z_index)

% Filter the data around the center frequency
tau = 0.2;
filter = exp(-tau*((Kx-x0).^2+(Ky-y0).^2+(Kz-z0).^2)); % the Gaussian function
path = zeros(3,49);
for i = 1:49
    rawData = reshape(subdata(:,i),n,n,n); 
    utn = fftshift(fftn(rawData));  
    unft = filter.*utn; % filter the data with the defined Gaussian function
    unf = ifftn(fftshift(unft)); % transform back to the spatial domain
    maxf = max(abs(unf),[],'all'); % Determine the locations of maximum frequency (the submarine)
    [xPath,yPath,zPath] = ind2sub(size(unf), find(abs(unf) == maxf));
    path(:,i) = [X(xPath,yPath,zPath), Y(xPath,yPath,zPath), Z(xPath,yPath,zPath)];
end

% Plot the trajectory of the submarine
plot3(path(1,:),path(2,:),path(3,:),'-o','Linewidth',1); grid on
set(gca, 'FontSize', 12)
title('Path of the submarine over time')
xlabel('x');
ylabel('y');
zlabel('depth');
axis([-10 10 -10 10 -10 10])
path(:,end) % the final location of the submarine