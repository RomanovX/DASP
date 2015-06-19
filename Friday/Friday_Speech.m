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
flag_bartp = true;
flag_sound = true;


%% Constants

SpT = 512;                                                      % Samples per time-frame
alpha = 0.95;                                                   % DD-weighting
beta = 0.8;                                                     % Bias compensation smoothing factor

overlap = 0.5;

var_est = 1;                                                    % Initial variance estimate (assuming white noise)
sp_est = 0;                                                     % Initial estimate of the clean speech signal

%% Loading audio & combining clean + noise

Fs = 16000;                                                     % Sampling frequency of 16000 Hz (same for all)   

[clean] = audioread('/audio/clean.wav');                        % Clean signal
[noise] = audioread('/audio/noise1.wav');                       % Noise

cleanpad = [clean; zeros(length(noise)-length(clean),1)];       % Zero padding

noisy = cleanpad + noise;                                       % Noisy 2

%% Signal Properties

    % Noise 1
    L_cleanpad = length(cleanpad);
    t_cleanpad = L_cleanpad/Fs;
    T_cleanpad = linspace(0,t_cleanpad, L_cleanpad)';

    L_noise = length(noise);
    t_noise = L_noise/Fs;
    T_noise = linspace(0,t_noise, L_noise)';

    L_noisy = length(noisy);
    t_noisy = L_noisy/Fs;
    T_noisy = linspace(0,t_noisy, L_noisy)';

if(flag_plots)  
    hold on
    subplot(4,1,1)
    plot(T_cleanpad, cleanpad)
    subplot(4,1,2)
    plot(T_noise, noise)
    subplot(4,1,3)
    plot(T_noisy, noisy)
end

%% Segmentation

Window = Modhanning(SpT);                                       % Modified Hanning curve
OS = fix(SpT*overlap);                                          % Calculate number of overlapping samples
N = fix((L_noisy-SpT)/OS +1);                                   % Number of segments

Index = (repmat(1:SpT,N,1)+repmat((0:(N-1))'*OS,1,SpT))';       % Index of overlapping samples
HW = repmat(Window,1,N);                                        % Hanning function copied for all segments
noisy_seg = noisy(Index).*HW;                                         % Apply Hanning function on each segment


%% DFT

Y = fft(noisy_seg);

% Note: this does a fft for every column, so for every time frame


%% Algorithm

% Bartlett Estimate

P_Y = (abs(Y).^2)/SpT;                                           % Periodogram per segment

% Computing the Bartlett estimate of the signal
Bart_Y = sum(P_Y,2) / N;

if(flag_bartp)
    figure
    plot(Bart_Y)
end
%%
% Paper method
h = waitbar(0, 'Waiting...');



for i = 2:N
   
    
end

waitbar(i/SpT,h,num2str(fix(100*i/SpT)));


for k=1:SpT
    
    var2(k,1) = var_est;
    S_est(k,1) = sp_est;
    
    for i=2:f
        magsY(k,i) = (abs(Y(k,i)))^2;
        SNR_ML(k,i) = max(((magsY(k,i)/(var2(k,i-1))) - 1), 0);             % Estimating SNR using ML
        aPost(k,i) = magsY(k,i)/var2(k,i-1);                                  % A Posteriori SNR
        PSD1(k,i) = ((1/((1+SNR_ML(k,i))^2)) + (SNR_ML(k,i)/((1+SNR_ML(k,i))*aPost(k,i)))) * magsY(k,i);
        SNR_DD(k,i) = alpha * ((abs(S_est(k,i-1)))^2)/(var2(k,i-1)) + (1-alpha) * max(((magsY(k,i)/(var2(k,i-1))) - 1), 0);
        Binv(k,i) = (1 + SNR_DD(k,i)) * gammainc(2, (1/(1+SNR_DD(k,i)))) + exp(-1/(1+SNR_DD(k,i)));
        B(k,i) = Binv(k,i)^(-1);
        var1(k,i) = PSD1(k,i) * B(k,i);
        var2(k,i) = beta * var2(k,i-1) + (1-beta) * var1(k,i);
        PSD2(k,i) = var2(k,i);
        
% Combining

        Hr(k,i) = 1 - (PSD2(k,i) / Bart_Y(k));
        S_est(k,i) = Hr(k,i) * Y(k,i);
        
        
        
    end
    
    
end

%% Inverse FFT

s_est = ifft(S_est);

%% De-Segment:

s = s_est(:);

%% Add to plot

if(flag_plots)

    L_s = length(s);
    t_s = L_s/Fs;
    T_s = linspace(0,t_s, L_s)'*1000;


    subplot(4,1,4)
    plot(T_s, s)
end


%% Sound output
if(flag_sound)
    sound(s,Fs)
end

close(h)
