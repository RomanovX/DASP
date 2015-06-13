%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     %
%             IN4182:DASP             %
%      Amritpal, Remy, Yadnyesh       %
%                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
close all

%% 1: Loading audio & combining clean + noise

Fs = 16000;                                                     % Sampling frequency of 16000 Hz (same for all)   

[clean] = audioread('/audio/clean.wav');                        % Clean signal
[noise1] = audioread('/audio/noise1.wav');                      % Noise 1
[noise2] = audioread('/audio/intersection_soundjay.wav');       % Noise 2

noisy1 = clean + noise1;
%noisy2 = clean + noise2;