function [p1] = OffsetCorrector(p1,up250)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
ang1=angle(p1(1:8*2048)./repmat(up250,8,1));
fOffset1=mean(removeoutliers(gradient(unwrap(ang1(512:512+1024)))));
for i=1:2048*55.25
    p1(i)=p1(i)*exp(-1j*i*fOffset1);
end
ang1=angle(p1(1:1*2048)./repmat(up250,1,1));
p=polyfit([1:1025]',unwrap(ang1(512:512+1024)),1);
for i=1:2048*55.25
    p1(i)=p1(i)*exp(-1j*i*p(1));
end
ang1=angle(p1(1:1*2048)./repmat(up250,1,1));
for i=1:2048*55.25
    p1(i)=p1(i)*exp(-1j*mean(ang1(512:512+1024)));
end

end

