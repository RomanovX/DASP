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

flag_plots = false;



%% Constants

SpT = 400;                                                      % Samples per time-frame

%% Loading audio & combining clean + noise

Fs = 16000;                                                     % Sampling frequency of 16000 Hz (same for all)   

[clean] = audioread('/audio/clean.wav');                        % Clean signal
[noise1] = audioread('/audio/noise1.wav');                      % Noise 1
[noise2] = audioread('/audio/intersection_soundjay.wav');       % Noise 2

cleanpad = [clean; zeros(length(noise2)-length(clean),1)];      % Zero padding

noisy1 = clean + noise1;                                        % Noisy 1
noisy2 = cleanpad + noise2;                                     % Noisy 2

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

f1 = ceil(L_noisy1/SpT);                                        % Number of frames Y1
f2 = ceil(L_noisy2/SpT);                                        % Number of frames Y2

y1 = zeros(512, f1);                                            % Zero-padding
y2 = zeros(512, f2);                                            % Zero-padding

for i = 1:f1-1
    y1(1 : SpT, i) = noisy1(SpT*(i-1)+1 : SpT*i, 1);            % Segmenting Y1
end
y1(1 : L_noisy1-((f1-1)*SpT)+1, f1) = noisy1(((f1-1)*SpT) : L_noisy1, 1); 
                                                                % Adding last frame  
for i = 1:f2-1
    y2(1 : SpT, i) = noisy2(SpT*(i-1)+1 : SpT*i, 1);            % Segmenting Y1
end
y2(1 : L_noisy2-((f2-1)*SpT)+1, f2) = noisy2(((f2-1)*SpT) : L_noisy2, 1); 
                                                                % Adding last frame 

%% DFT

Y1 = fft(y1);
Y2 = fft(y2);

% Note: this does a fft for every column, so for every time frame

%% A priori SNR using ML



%% PSD Estimation

%% A priori SNR using DD

%% Bias compensation

%% Smoothing

%% Locking prevention ?????????

