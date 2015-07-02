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
flag_bartp = false;
flag_sound = true;


%% Constants

SpT = 512;                                                      % Samples per time-frame
alpha = 0.95;                                                   % DD-weighting
beta = 0.8;                                                     % Bias compensation smoothing factor

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
    T_cleanpad = linspace(0,t_cleanpad, L_cleanpad)'*1000;

    L_noise = length(noise);
    t_noise = L_noise/Fs;
    T_noise = linspace(0,t_noise, L_noise)'*1000;

    L_noisy = length(noisy);
    t_noisy = L_noisy/Fs;
    T_noisy = linspace(0,t_noisy, L_noisy)'*1000;

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

f = L_noisy/SpT;                                                % Number of frames Y

y = zeros(SpT,f);

for j = 1:f
    y(1 : SpT, j) = noisy(SpT*(j-1)+1 : SpT*j, 1);              % Segmenting Y
end



%% DFT

Y = fft(y);

theta=angle(Y);

% Note: this does a fft for every column, so for every time frame


%% Algorithm

% Bartlett Estimate

P_Y = (abs(Y).^2)/SpT;                                           % Periodogram per segment

% Computing the Bartlett estimate of the signal
Bart_Y = sum(P_Y,2) / f;

if(flag_bartp)
    figure
    plot(Bart_Y)
end

% Paper method
h = waitbar(0, 'Waiting...');

for k=1:SpT
    
%     var2(k,1) = var_est;
    S_est(k,1) = sp_est;
    PSD2(k,1) = var_est / SpT;
    
    for j=2:f
        magsY(k,j) = (abs(Y(k,j)))^2;
        SNR_ML(k,j) = max(((magsY(k,j)/(PSD2(k,j-1))) - 1), 0);             % Estimating SNR using ML [Hendriks MMSE]
        aPost(k,j) = magsY(k,j)/PSD2(k,j-1);                                  % A Posteriori SNR
        PSD1(k,j) = ((1/((1+SNR_ML(k,j))^2)) + (SNR_ML(k,j)/((1+SNR_ML(k,j))*aPost(k,j)))) * magsY(k,j);
        SNR_DD(k,j) = alpha * ((abs(S_est(k,j-1)))^2)/(PSD2(k,j-1)) + (1-alpha) * max(((magsY(k,j)/(PSD2(k,j-1))) - 1), 0);
        Binv(k,j) = (1 + SNR_DD(k,j)) * gammainc(2, (1/(1+SNR_DD(k,j)))) + exp(-1/(1+SNR_DD(k,j)));
%         B(k,j) = Binv(k,j)^(-1);
        PSD1(k,j) = PSD1(k,j)/Binv(k,j);
        PSD2(k,j) = beta * PSD2(k,j-1) + (1-beta) * PSD1(k,j);
%         PSD2(k,i) = var2(k,i);
        
% Combining

        Hr(k,j) = max((1 - (PSD2(k,j) / Bart_Y(k))),0);
        S_est(k,j) = Hr(k,j) *abs(Y(k,j));
        S_est(k,j)=S_est(k,j)*exp(i*theta(k,j));
        
        
    end
    
    waitbar(k/SpT,h,num2str(fix(100*k/SpT)));
end

%% Inverse FFT

s_est = real(ifft(S_est));

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
