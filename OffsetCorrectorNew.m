function [sig] = OffsetCorrectorNew(sig,ref)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%ang1=angle(p1(1:8*2048)./repmat(up250,8,1));
fOffset1=mean(removeoutliers(gradient(unwrap(angle(sig./ref)))));
for i=1:2048*55.25
    sig(i)=sig(i)*exp(-1j*i*fOffset1);
end
ang1=angle(sig(1:1*2048)./ref(1:1*2048));
p=polyfit([1:1025]',unwrap(ang1(512:512+1024)),1);
for i=1:2048*55.25
    sig(i)=sig(i)*exp(-1j*i*p(1));
end
p=polyfit([1:length(angle(sig./ref))]',unwrap((angle(sig./ref))),1);
for i=1:2048*55.25
    sig(i)=sig(i)*exp(-1j*i*p(1));
end
ang1=angle(sig./ref);
if(max(ang1(1:1e5))-min(ang1(1:1e5))>6)
    ang1=mod(ang1,2*pi);
   for i=1:2048*55.25
        sig(i)=sig(i)*exp(-1j*mean(ang1(1:1e5)));
    end
else
    for i=1:2048*55.25
        sig(i)=sig(i)*exp(-1j*mean(ang1(1:1e5)));
    end
    
end

end

