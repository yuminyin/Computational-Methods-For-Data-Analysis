%% guitar in the GNR clip
clear all; close all; clc

[y, Fs] = audioread('GNR.m4a');
tr_gnr = length(y)/Fs; % record time in seconds

S = y'; % Signal
L = tr_gnr;
n = length(S);
t2 = linspace(0,L,n+1); t = t2(1:n); % vector of time points
k = (1/L)*[0:n/2-1 -n/2:-1]; ks = fftshift(k); % Notice the 1/L instead of 2*pi/L

% Apply the Gabor transform implementation
tau = 0:0.1:L; %the center of the window
a = 1000;
for i = 1:length(tau)
    g = exp(-a*(t-tau(i)).^2); % Window function
    Sg = g.*S;
    Sgt = fft(Sg);
    Sgt_spec(:,i) = fftshift(abs(Sgt));
end

% construct the spectrum
pcolor(tau,ks,Sgt_spec) 
shading interp
set(gca, 'ylim', [0 1000], 'FontSize', 14)
xlabel('time (t)'), ylabel('frequency (k)')
colormap(hot)
yticks([0 278 311 370 415 500 554 698 1000]);
yticklabels({'0', 'C^#', 'D^#', 'F^#', 'G^#', '500', 'C^#', 'F', '1000'});
title('Spectrum of the guitar in the GNR clip', 'Fontsize', 14)

%% isolate the bass in Comfortably Numb
clear all; close all; clc

[y, Fs] = audioread('Floyd.m4a');
tr_gnr = length(y)/Fs; % record time in seconds

S = y'; % Signal
L = tr_gnr;
n = length(S);
t2 = linspace(0,L,n+1); t = t2(1:n); % vector of time points
k = (1/L)*[0:n/2-1 -n/2:-1]; ks = fftshift(k); % Notice the 1/L instead of 2*pi/L

% bandpass filter
S_fft = fft(S);
S_filter = S_fft.*fftshift(60<abs(S_fft)<250);
S_bass = ifft(S_filter);

% Apply the Gabor transform implementation
tau = 0:2.5:L; %the center of the window
a = 10;
for i = 1:length(tau)
    g = exp(-a*(t-tau(i)).^2); % Window function
    Sg = g.*S_bass;
    Sgt = fft(Sg);
    Sgt_spec(:,i) = fftshift(abs(Sgt));
end

% construct the spectrum
Sgt_spec(end,:) = [];
pcolor(tau,ks,Sgt_spec) 
shading interp
set(gca, 'ylim', [60 250], 'FontSize', 16)
xlabel('time (t)'), ylabel('frequency (k)')
colormap(hot)
title('Spectrum of the bass in the Floyd clip', 'Fontsize', 16)

%% isolate 10s of the bass in Comfortably Numb
clear all; close all; clc

[y, Fs] = audioread('10sFloyd.m4a');
tr_gnr = length(y)/Fs; % record time in seconds

S = y'; % Signal
L = tr_gnr;
n = length(S);
t2 = linspace(0,L,n+1); t = t2(1:n); % vector of time points
k = (1/L)*[0:n/2-1 -n/2:-1]; ks = fftshift(k); % Notice the 1/L instead of 2*pi/L

% bandpass filter
S_fft = fft(S);
S_filter = S_fft.*fftshift(60<abs(S_fft)<250);
S_bass = ifft(S_filter);

% Apply the Gabor transform implementation
tau = 0:0.25:L; %the center of the window
a = 10;
for i = 1:length(tau)
    g = exp(-a*(t-tau(i)).^2); % Window function
    Sg = g.*S_bass;
    Sgt = fft(Sg);
    Sgt_spec(:,i) = fftshift(abs(Sgt));
end

% construct the spectrum
Sgt_spec(end,:) = [];
pcolor(tau,ks,Sgt_spec) 
shading interp
set(gca, 'ylim', [60 250], 'FontSize', 16)
xlabel('time (t)'), ylabel('frequency (k)')
colormap(hot)
yticks([60 81 89 97 110 123 250]);
yticklabels({'60', 'E', 'F', 'G', 'A', 'B', '250'});
title('Spectrum of the 10s bass in the Floyd clip', 'Fontsize', 16)

%% the Guitar in Comfortably Numb
clear all; close all; clc

[y, Fs] = audioread('Floyd.m4a');
tr_gnr = length(y)/Fs; % record time in seconds

S = y'; % Signal
L = tr_gnr;
n = length(S);
t2 = linspace(0,L,n+1); t = t2(1:n); % vector of time points
k = (1/L)*[0:n/2-1 -n/2:-1]; ks = fftshift(k); % Notice the 1/L instead of 2*pi/L

% bandpass filter
S_fft = fft(S);
S_filter = S_fft.*fftshift(200<abs(S_fft)<500);
S_guitar = ifft(S_filter);

% Apply the Gabor transform implementation
tau = 0:2.5:L; %the center of the window
a = 10;
for i = 1:length(tau)
    g = exp(-a*(t-tau(i)).^2); % Window function
    Sg = g.*S_guitar;
    Sgt = fft(Sg);
    Sgt_spec(:,i) = fftshift(abs(Sgt));
end

% construct the spectrum
Sgt_spec(end,:) = [];
pcolor(tau,ks,Sgt_spec) 
shading interp
set(gca, 'ylim', [200 500], 'FontSize', 16)
xlabel('time (t)'), ylabel('frequency (k)')
colormap(hot)
title('Spectrum of the guitar in the Floyd clip', 'Fontsize', 16)

%% the 10s of Guitar in Comfortably Numb
clear all; close all; clc

[y, Fs] = audioread('10sFloyd.m4a');
tr_gnr = length(y)/Fs; % record time in seconds

S = y'; % Signal
L = tr_gnr;
n = length(S);
t2 = linspace(0,L,n+1); t = t2(1:n); % vector of time points
k = (1/L)*[0:n/2-1 -n/2:-1]; ks = fftshift(k); % Notice the 1/L instead of 2*pi/L

% bandpass filter
S_fft = fft(S);
S_filter = S_fft.*fftshift(200<abs(S_fft)<500);
S_guitar = ifft(S_filter);

% Apply the Gabor transform implementation
tau = 0:0.5:L; %the center of the window
a = 10;
for i = 1:length(tau)
    g = exp(-a*(t-tau(i)).^2); % Window function
    Sg = g.*S_guitar;
    Sgt = fft(Sg);
    Sgt_spec(:,i) = fftshift(abs(Sgt));
end

% construct the spectrum
Sgt_spec(end,:) = [];
pcolor(tau,ks,Sgt_spec) 
shading interp
set(gca, 'ylim', [200 500], 'FontSize', 16)
xlabel('time (t)'), ylabel('frequency (k)')
colormap(hot)
yticks([200 220 247 261 294 330 370 392 440 500]);
yticklabels({'200', 'A', 'B', 'C', 'D', 'E', 'F^#', 'G', 'A', '500'});
title('Spectrum of the guitar in the Floyd clip', 'Fontsize', 16)
