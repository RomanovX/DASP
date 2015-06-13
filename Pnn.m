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

%% DFT

%% A priori SNR using ML

%% PSD Estimation

%% A priori SNR using DD

%% Bias compensation

%% Smoothing

%% Locking prevention ?????????

