function out = simple_bandpass(in, pct, pct2)
% This is a sharp bandpass that only keeps signals between pct% to pct2% of
% the entire spectrum
    n = length(in);
    a = fft(in);
    
    % -f ____ -cutoff ^^^^ -cutoff2 ___________ cutoff2 ^^^^ cutoff ___ f
    cutoff = round((pct/100)*n/2);
    cutoff2 = round((pct2/100)*n/2);
    midpoint = round(n/2);
%     disp("bandpass params")
%     disp([midpoint, cutoff, cutoff2]);
%     disp([1 ,cutoff, midpoint-cutoff2, midpoint+cutoff2, n-cutoff, n])
    a(1:cutoff) = 0;
    a(end-cutoff:end)=0;
    a(midpoint-cutoff2:midpoint+cutoff2)=0;
    out = ifft(a);
end