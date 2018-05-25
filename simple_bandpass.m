function out = simple_bandpass(in, startPercent, endPercent)
	%%% simple bandpass removes out all frequency components outside
	%%% startPercent and endPercent  
	n = length(in);
	a = fft(in);

	cutoff = round((startPercent/100)*n/2);
	cutoff2 = round((endPercent/100)*n/2);
	midpoint = round(n/2);
	
	a(1:cutoff) = 0;
	a(end-cutoff:end)=0;
	a(midpoint-cutoff2:midpoint+cutoff2)=0;
	
	out = ifft(a);
end