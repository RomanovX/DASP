function [psd_plot] = plot_psd (psd, fs, framelength)
% Plots PSD spectrogram in dB, in the latest opened figure() window.
% Also returns the datamatrix that was actually plotted (fft-shifted and dB-converted data).
% You need to insert a title and colorbar yourself.
% Input arguments:
% fs is the used sampling rate, framelength is the framelength used
% for framing the the time signal (NOT THE SAME AS NFFT!).

	% Find parameters
	[NFFT, frames] = size(psd);

	% Define frequency axis
	f = linspace(-fs/2, fs/2, NFFT);

	% FFTshift and convert to dB
	psd = fftshift (abs(psd), 1);
	psd_plot = 10*log10(psd);

	% Plot results
	imagesc(1:frames, f, psd_plot);
	xlabel('Frame index');
	ylabel('Frequency bin [Hz]');
	zlabel('Magnitude [dB]');
%	caxis manual;
%	caxis([minimum maximum]);
	colorbar;
