function window = mod_hanning(length);

    error(nargchk(1,1,nargin));
    window =(0.5 + 0.5*cos(2*pi*(-(length-1)/2:(length-1)/2)/length))';

end