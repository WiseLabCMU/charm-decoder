function [ SR ] = calcSR( packet,down )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    maxe=0;
    for i=1:length(packet)-length(down)-1
        t=max(abs(fftshift(fft(packet(i:i+length(down)-1).*down))));
        if(t>maxe)
            maxe=t;
        end
    end
    SR=maxe^2;
end

