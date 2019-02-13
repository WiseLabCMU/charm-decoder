function out = fftBandpass(x, Wstart, Wend)
%FFTBANDPASS A sharp Bandpass filter that filters out everything outside of
% WStart and Wend
%   Detailed explanation goes here

    n = length(in);
    a = fft(in);
    
    % -f ____ -cutoff ^^^^ -cutoff2 ___________ cutoff2 ^^^^ cutoff ___ f
    midpoint = round(n/2);
    cutoffLow = round(Wstart*n/2);
    cutoffHigh = round(Wend*n/2);
    
    a(1:cutoffLow) = 0;
    a(midpoint-cutoffHigh:midpoint+cutoffHigh)=0;
    a(end-cutoffLow:end)=0;
    out = ifft(a);
end

