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

IS = 1;                                                         % Seconds of initial silence (no present speech)



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
noisy_seg = noisy(Index).*HW;                                   % Apply Hanning function on each segment


%% DFT

Y = fft(noisy_seg);
% Note: this does a fft for every column, so for every time frame

% Spectogram:

%Y_phase257 = angle((Y(1:fix(end/2)+1,:)));                       % Only looking at half because after is complex conjugate            
Y_phase512 = angle(Y);
%Y_real257 = abs(Y(1:fix(end/2)+1,:)); 
Y_real512 = abs(Y);

magsY = Y_real512.^2;                                            % Square of magnitude

%% Initializing

IS_seg = fix((IS*Fs-SpT)/OS +1);                                 % Segments of Initial Silence (same formula as in Segmentation #3                                           

NPS_mean = mean(Y_real512(:,1:IS_seg),2);                        % Mean of coeffs during speech absence
NPS_var = var(Y_real512(:,1:IS_seg),0,2);                        % Variance of coeffs during speech absence

PSD = zeros(size(Y));                                            % Initializing PSD
PSD(:,1) = (NPS_var/SpT);                                        % Using variance/SpT as initial estimate

aPost = ones(SpT,1);                                             % Initial value for a posteriori SNR

S = zeros(size(Y));                                              % Assuming there is no speech in the initial segment


%% Bartlett Estimate

% P_Y = (abs(Y).^2)/SpT;                                           % Periodogram per segment
% 
% % Computing the Bartlett estimate of the signal
% Bart_Y = sum(P_Y,2) / N;
% 
% if(flag_bartp)
%     figure
%     plot(Bart_Y)
% end

Welch = pwelch(noisy_seg,SpT,OS, SpT);
Welch_flip = flipud(Welch);
Welch512 = [Welch(1:end-1,:); Welch_flip(2:end, :)];

Welch_avg = mean(Welch512,2);


%% Paper method
h = waitbar(0, 'Waiting...');

for i = 2:N
    varN = PSD(:,i-1)*SpT;
    SNR_ML = max(((magsY(:,i)./(PSD(:,i-1))) - 1), 0);                                   % Estimating a priori SNR using ML
    aPost_new = magsY(:,i)./varN;                                                     % New a posteriori SNR (Assuming Noise PSD is relatively constant)
    PSD_MMSE = ((1./((1+SNR_ML).^2)) + (SNR_ML./((1+SNR_ML).*aPost))) .* magsY(:,i);
    SNR_DD = alpha*(abs(S(:,i-1)).^2./varN)+(1-alpha).*max(aPost_new-1,0);                          %Decision Directed Method for A Priori SNR
    aPost = aPost_new;                                                 % Retains a posteriori SNR of current time frame for next
    B = (1 + SNR_DD) .* gammainc(2, (1./SNR_DD)) + exp(-1./(1+SNR_DD));
    B = B.^(-1);
    BiasComp = PSD_MMSE .* B;
    PSD(:,i) = beta * PSD(:,i-1) + (1-beta) * BiasComp;
    
    H = 1 - PSD(:,i)./Welch_avg;
    
    S(:,i) = H.*Y(:,i);
    %S(:,i) = S(:,i).*exp(1j*Y_phase512(:,i));    
  
    
    waitbar(i/N,h,num2str(fix(100*i/N)));
    
end

close(h);

%% De-segment and IFFT

stemp = ifft(S);

%Spec=S.*exp(1i * Y_phase512);

%Spec=S.*exp(1i*Y_phase512);

%if mod(SpT,2) %if FreqResol is odd
%     Spec=[Spec;flipud(conj(Spec(2:end,:)))];
% else
%     Spec=[Spec;flipud(conj(Spec(2:end-1,:)))];
% end

s = zeros((N-1)*OS+SpT,1);

for i=1:N
    start = (i-1)*OS+1;    
    %s(start:start+SpT-1)=s(start:start+SpT-1)+real(ifft(Spec(:,i),SpT));  
    s(start:start+SpT-1)=s(start:start+SpT-1)+stemp(:,i); 
end


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

%close(h)
