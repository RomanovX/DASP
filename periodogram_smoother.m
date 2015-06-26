
% X: framed time data
function Px = periodogram_smoother (X, NFFT)
% Smoothed periodogram over the windowed frames in X, using
% an exponential recursive smoothing window.
% Each column in X should represent a time-frame.
% ASSUMES THAT MOD_HANNING IS USED FOR WINDOWING!

	alpha = 0.9;	% Smoothing parameter

	prev_psd = zeros(NFFT,1);

	[framelength, nof_frames] = size(X);

	window = mod_hanning(framelength);
	U = 1/framelength * sum(window.^2); 

	for frame=1:nof_frames
		x = X(:,frame);
		Px_frame = alpha*prev_psd + (1-alpha)*abs(fft(x, NFFT)).^2;
		prev_psd = Px_frame;
		Px(:,frame) = Px_frame;
	end;
	Px = (1/U)*(1/framelength)*Px;
    
end

