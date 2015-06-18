%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     %
%             IN4182:DASP             %
%      Amritpal, Remy, Yadnyesh       %
%                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Clearing
    clear all
    clc
    close all
    
%% Flags

flag_plots = true;



%% Constants

SpT = 512;                                                      % Samples per time-frame
alplha = 0.98;                                                  % DD-weighting

%% Loading audio & combining clean + noise

Fs = 16000;                                                     % Sampling frequency of 16000 Hz (same for all)   

[clean] = audioread('/audio/clean.wav');                        % Clean signal
[noise1] = audioread('/audio/noise1.wav');                      % Noise 1
[noise2] = audioread('/audio/intersection_soundjay.wav');       % Noise 2

cleanpad = [clean; zeros(length(noise2)-length(clean),1)];      % Zero padding

noisy1 = clean + noise1;                                        % Noisy 1
noisy2 = cleanpad + noise2;                                     % Noisy 2

noisy1 = [noisy1; zeros(SpT-mod(length(noisy1),SpT),1)];        % Zero pad to length, multiple of 512
noisy2 = [noisy2; zeros(SpT-mod(length(noisy2),SpT),1)];        % Zero pad to length, multiple of 512

%% Signal Properties

    % Noise 1
    L_clean = length(clean);
    t_clean = L_clean/Fs;
    T_clean = linspace(0,t_clean, L_clean)'*1000;

    L_noise1 = length(noise1);
    t_noise1 = L_noise1/Fs;
    T_noise1 = linspace(0,t_noise1, L_noise1)'*1000;

    L_noisy1 = length(noisy1);
    t_noisy1 = L_noisy1/Fs;
    T_noisy1 = linspace(0,t_noisy1, L_noisy1)'*1000;

    % Noise 2
    L_cleanpad = length(cleanpad);
    t_cleanpad = L_cleanpad/Fs;
    T_cleanpad = linspace(0,t_cleanpad, L_cleanpad)'*1000;

    L_noise2 = length(noise2);
    t_noise2 = L_noise2/Fs;
    T_noise2 = linspace(0,t_noise2, L_noise2)'*1000;

    L_noisy2 = length(noisy2);
    t_noisy2 = L_noisy2/Fs;
    T_noisy2 = linspace(0,t_noisy2, L_noisy2)'*1000;

if(flag_plots)    
    figure
    subplot(3,2,1)
    plot(T_clean, clean)
    subplot(3,2,3)
    plot(T_noise1, noise1)
    subplot(3,2,5)
    plot(T_noisy1, noisy1)
    
    subplot(3,2,2)
    plot(T_cleanpad, cleanpad)
    subplot(3,2,4)
    plot(T_noise2, noise2)
    subplot(3,2,6)
    plot(T_noisy2, noisy2)
end

%% Segmentation

f1 = L_noisy1/SpT;                                              % Number of frames Y1
f2 = L_noisy2/SpT;                                              % Number of frames Y2

for i = 1:f1-1
    y1(1 : SpT, i) = noisy1(SpT*(i-1)+1 : SpT*i, 1);            % Segmenting Y1
end
for i = 1:f2-1
    y2(1 : SpT, i) = noisy2(SpT*(i-1)+1 : SpT*i, 1);            % Segmenting Y2
end

%% DFT

Y1 = fft(y1);
Y2 = fft(y2);

% Note: this does a fft for every column, so for every time frame

%% Paper method


for k=1:SpT
    for i=1:f1
SNR1(



% Starting SNR value
SNR1 = zeros(DFT,1);                                         
SNR2 = zeros(DFT,1);

Est1 = zeros(DFT,1);
Est2 = zeros(DFT,1);

Y1mag = abs(Y1).^2;
Y2mag = abs(Y2).^2;

var1 = var(Y1);
var2 = var(Y2);

apost1 = Y1mag/var1;

PrevNPSD = ???
% Estimator

for i = 2:f1
    SNR_ML(i) = max((Y1mag(i) / prev
    Est1(i) = ( (1/(1 + SNR(i))) + SNR(i)/((1 + SNR(i))*apost1;
    SNR_DD(i) = alpha*??? + (1-alpha)*max(Y1mag(i-1) - PrevNPSD, 0)
end








%% PSD Estimation

%% A priori SNR using DD

%% Bias compensation

%% Smoothing

%% Locking prevention ?????????

