%% QUESTION 1
% DO NOT REMOVE THE LINE BELOW
% MAKE SURE 'eel3135_lab11_comment.m' IS IN THE SAME DIRECTORY AS THIS FILE
clear; close all; clc;
type('eel3135_lab10_comment.m')

%% Question 2
n=0:59;
x = 0.75 + cos(pi*n/20) + cos(pi*n/15) + cos(pi*n + 2*pi/3);
w_DTFT = linspace(0, 2*pi-pi/5000, 10000);

% NOTE: USE THE FOLLOWING COMMENTED LINE FOR PLOTTING THE DFT ATOP THE DTFT
% (YOU NEED TO DEFINE w_DFT), THIS WILL MAKE THE PLOTS EASIER TO INTERPRET
% plot(w_DTFT,abs(X_DTFT)); 
% hold on; plot(w_DFT,abs(X_DFT),'.', 'markersize', 10); hold off;

%% Question 2a - Changes in the code below

%% Question 2b

x_dft = DFT(x);
f = (0:length(x_dft)-1)*(1/length(x_dft)); % frequency bins(normalized)
figure;
plot(f, abs(x_dft));
xlabel('Normalized Frequency (cycles/sample)');
ylabel('Magnitude');
title('Magnitude of DFT of x[n]');
grid on;

x_dtft = DTFT(x, w_DTFT);
figure;
plot(w_DTFT, abs(x_dtft));
xlabel('Normalized Frequency \omega (Radians/Sample)');
ylabel('Magnitude');
title('Magnitude of DTFT of x[n] for 60 samples');
grid on;

%% Question 2c

% Frequency vector for DFT (normalized)
w_dft = 2*pi*(0:length(x)-1)/length(x);

figure;
x_dtft = DTFT(x, w_DTFT);
plot(w_DTFT, abs(x_dtft));
x_dft = DFT(x);
hold on;
plot(w_dft, abs(x_dft),'.', 'MarkerSize',10); % Use red circles for DFT points
legend('DTFT', '60-length DFT');
hold off;
xlabel('Normalized Freqeuncy \omega (rads/s)');
ylabel('Magnitude');
title('DFT and DTFT Plot for 60 samples');
grid on;

%% Question 2d

figure;
x_1 = x(1:55);
% Freqency vector for DFT (normalized)
w_dft = 2*pi*(0:length(x_1)-1)/length(x_1);
x_dtft = DTFT(x_1, w_DTFT);
plot(w_DTFT, abs(x_dtft));
x_dft = DFT(x_1);
hold on;
plot(w_dft, abs(x_dft),'.', 'MarkerSize',10); % Use red circles for DFT points
legend('DTFT', '55-length DFT');
hold off;
xlabel('Normalized Freqeuncy \omega (rads/s)');
ylabel('Magnitude');
title('DFT and DTFT Plot for 55 samples');
grid on;

%% Question 2e

figure;
n = 0:64;
x = 0.75 + cos(pi*n/20) + cos(pi*n/15) + cos(pi*n + 2*pi/3);
% Freqency vector for DFT (normalized)
w_dft = 2*pi*(0:length(x)-1)/length(x);
x_dtft = DTFT(x, w_DTFT);
plot(w_DTFT, abs(x_dtft));
x_dft = DFT(x);
hold on;
plot(w_dft, abs(x_dft),'.', 'MarkerSize',10); % Use red circles for DFT points
legend('DTFT', '65-length DFT');
hold off;
xlabel('Normalized Freqeuncy \omega (rads/s)');
ylabel('Magnitude');
title('DFT and DTFT Plot for 65 samples');
grid on;

%% Question 2f

figure;
n = 0:199;
x = 0.75 + cos(pi*n/20) + cos(pi*n/15) + cos(pi*n + 2*pi/3);
% Freqency vector for DFT (normalized)
w_dft = 2*pi*(0:length(x)-1)/length(x);
x_dtft = DTFT(x, w_DTFT);
plot(w_DTFT, abs(x_dtft));
x_dft = DFT(x);
hold on;
plot(w_dft, abs(x_dft),'.', 'MarkerSize',10); % Use red circles for DFT points
legend('DTFT', '200-length DFT');
hold off;
xlabel('Normalized Freqeuncy \omega (rads/s)');
ylabel('Magnitude');
title('DFT and DTFT Plot for 200 samples');
grid on;

%% Question 2g

% The DF is the same as the DFT, but it's in discrete time. The DFT is
% sampling a few points from the DTFT. This makes sense as instead of a
% continuous omega_hat term, the DFT transform uses the discreter term
% 2*pi*k/N. It is easier for computers to utilize since it's sampling
% the DTFT with N samples. The result would be the same when there is
% double over the number of samples as the number of input values.

%% Question 3
clear all; clc;
x = [zeros(1,10) ones(1,35)];


%% Question 3a - At the function of for IDFT

%% Question 3b

n = 0:99;
x_2 = [x zeros(1, 55)];
figure;
stem(n, x_2);
xlabel('Samples');
ylabel('x[n]');
title('Stem Plot of Input Signal');

%% Question 3c

y_n = conv(x_2, x_2);
y_sampled = y_n(1:100);
figure;
stem(n, y_sampled);
xlabel('Samples');
ylabel('y[n]');
title('Stem Plot of Output Signal');

%% Question 3d

x_dft_2 = DFT(x_2); % Converts to time domain
% Multiply both DFT's together to get frequency domain equivalent of
% convolution
x_s = x_dft_2.*x_dft_2;

% Take inverse DFT
y = IDFT(x_s);

figure;
stem(n, y);
xlabel('Samples');
ylabel('y[n]');
title('Stem Plot of Output Signal using DFT and IDFT');
ylim([0,35]);

%% Question 3e

n2 = 0:59;
x_60 = [x zeros(1,15)];
x_60 = DFT(x_60);
x_s60 = x_60.*x_60;
y_60 = IDFT(x_s60);

figure;
stem(n2, y_60(1:60));
xlabel('Samples');
ylabel('y[n]');
title('Stem Plot of Output Signal using DFT and IDFT for 60 Samples');


%% Question 3f

% There is a big difference between the last two solutions and this is due
% to the zero-padding. This could be because DFT and IDFT rely on circular
% convolution, is the length of the sampled signal is shorter than the
% length of the of the convolution, then some previous valies would wrap
% around and affect DFT (ailiasing). This happens with the length 60 input
% signal since it is not 2*length(x). However, the 100 length signal avoids
% this because it is padded to a value greater than 90, so it never wraps
% around. 

%% Question 4
% Choose a song at least 3 minutes long to use in this
% problem, and include it in your submission. Load it into MATLAB using
% audioread. Note that mose audio files will be stereo, so you need to make
% Sure that you only use one column of audio data for this part of the lab
% This site has a large archive of free music that you can choose from:
% https://freemusicarchive.org/static
clear all; clc;
[x,fs] = audioread('porcelain-face-official-audio.wav');
x = x(:,1);

%% Question 4a

x3 = x(1:10000);
n = 0:9999;
w_DFT = 2*pi*(0:length(x3)-1)/length(x3);

tic;
x3 = DFT(x3);
time = toc;

disp(['The calculation time is: ', num2str(time), ' seconds']);
figure;
plot(w_DFT, abs(x3));
xlabel('Normalized Freqeuncy \omega (rads/s)');
ylabel('Magnitude');
title('DFT Plot for Audio Sample');


%% Question 4b

x3_f = x(1:10000);
n = 0:9999;
w_DFT = 2*pi*(0:length(x3_f)-1)/length(x3_f);

tic;
x3_f = fft(x3_f, length(x3_f));
time = toc;

disp(['The calculation time is: ', num2str(time), ' seconds']);
figure;
plot(w_DFT, abs(x3_f));
xlabel('Normalized Freqeuncy \omega (rads/s)');
ylabel('Magnitude');
title('FFT Plot for Audio Sample');

%% Question 4c

% The difference in the magnitude results aren't apparent. The audio sample
% is relatively long since it has 10,000 samples, but both magnitude plots
% appear the same even though they were calculated differently. 

%% Question 4d

% THe FFT is much faster than the DFT. 

%% Question 4e

w_DFT = 2*pi*(0:length(x)-1)/length(x);

tic;
x = fft(x);
time = toc;

disp(['The calculation time is: ', num2str(time), ' seconds']);
figure;
plot(w_DFT, abs(x));
xlabel('Normalized Freqeuncy \omega (rads/s)');
ylabel('Magnitude');
title('FFT Plot for Audio Sample');

%% Functions provided for the lab
function H = DTFT(x,w)
% DTFT(X,W)  compute the Discrete-time Fourier Transform of signal X
% acroess frequencies defined by W. 

    H = zeros(1, length(w));
    for nn = 1:length(x)
        H = H + x(nn).*exp(-1j*w.*(nn-1));
    end
    
end

function X = DFT(x)
% DFT(x)  compute the N-point Discrete Fourier Transform of signal x  
% Where N is the length of signal x
    N = length(x); % FILL THIS LINE IN
    w = 2*pi*(0:N-1)/N; % FILL THIS LINE IN
    X = zeros(1, length(w));    
    for nn = 1:length(x)    
        X = X + x(nn).*exp(-1j*w.*(nn-1));    
    end
end

function x = IDFT(X)
% IDFT(x)  compute the N-point Inverse Discrete Fourier Transform of signal
% X where N is the length of signal X
    N = length(X);% FILL THIS LINE IN
    w = 2*pi*(0:N-1)/N;% FILL THIS LINE IN
    x = zeros(1, length(w));    
    for nn = 1:length(X)    
        x = x + X(nn)*exp(1j*w.*(nn-1));% FILL THIS LINE IN    
    end
    x=x/N;
end