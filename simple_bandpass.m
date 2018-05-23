function out = simple_bandpass(in, pct, pct2)
    n = length(in);
    a = fft(in);
    cutoff = round((pct/100)*n/2);
    cutoff2 = round((pct2/100)*n/2);
    midpoint = round(n/2);
    a(1:cutoff) = 0;
    a(end-cutoff:end)=0;
    a(midpoint-cutoff2:midpoint+cutoff2)=0;
    out = ifft(a);
end