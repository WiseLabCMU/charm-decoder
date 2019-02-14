function p1 = getPacket(sig1,down250,up250)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
[y,x]=xcorr(sig1,repmat(down250,2,1));
plot(x,abs(y));
[~,locs]=findpeaks(abs(y),'MinPeakProminence',0.3,'sortstr','descend');
loc1=x(median(locs(1:3)));
MASTER=[];
for i=loc1-2048*10:loc1-2048*10+100
    [y,x]=xcorr(sig1(i:i+2048),up250);
    MASTER=[MASTER max(abs(y))];
end
l = find(MASTER(2:101)-MASTER(1:100) < 0, 1);
p1=sig1(loc1-2048*10+l:loc1+2048*45.25+l);

end

