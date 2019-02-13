tic

BW = 1.25e5; % in Hz
% SF = 10; % spreading factor: 
% NOTE:Unused
chirp_size = 576; % chirp length in number of samples
% TODO: change this to a known constant
packet_length = 20; % packet size in number of chirps
packet_size = chirp_size*packet_length;
bandwidth_sampling_factor = 1; % table  Sample/KHz       bandwidth_sampling_factor
%         125                1
%         250                2
%         500                3
Fs = BW;
symbol_length = chirp_size;
symbol_length_upsampled = bandwidth_sampling_factor*chirp_size;
freq_shift_per_sample =  Fs/symbol_length; % How each frequency bin maps to a difference in frequency
Ts = 1/freq_shift_per_sample; % Symbol Duration
f = linspace(-BW/2, BW/2-freq_shift_per_sample, symbol_length); % The X-Axis
% of the FFT plot
reset_freq = -BW/2; % The initial frequency of the base chirp
final_freq = (BW/2);%-freq_shift_per_sample; % The final frequency
[up,down] = my_create_chirpspecial1(bandwidth_sampling_factor*Fs,Ts,reset_freq,final_freq,chirp_size);

disp("created reference chirps")
toc

nsym=38.25;
len = zeros(1, 8);
rdata_struct = cell(1, 8);
for ii=1:8
    rdata_struct{ii}=importdata(sprintf('data/p1_sf7_%d_IQ.csv', ii));
    len(ii) = length(rdata_struct{ii});
end

disp("loaded data")
toc

ps_struct = [];

data=rdata_struct{5};
    temp=1i*data(:, 1)+data(:, 2);
    temp=simple_bandpass(temp, 2, 0);
    maxi=0;
    maxval=0;
    for i=packet_size+1:1000:length(temp)-packet_size-1
        if(sumsqr(real(temp(i:i+packet_size)))/sumsqr(real(temp(i-packet_size:i)))>maxval)
            maxi=i;
            maxval=sumsqr(real(temp(i:i+packet_size)))/sumsqr(real(temp(i-packet_size:i)));
        end
    end
    for i=maxi-1000:1:maxi-1000
        if(sumsqr(real(temp(i:i+packet_size)))/sumsqr(real(temp(i-packet_size:i)))>maxval)
            maxi=i;
            maxval=sumsqr(real(temp(i:i+packet_size)))/sumsqr(real(temp(i-packet_size:i)));
        end
    end
    rval2=maxi;

loc_weight=[];
data=rdata_struct{1};
    temp=1i*data(:, 1)+data(:, 2);
    temp=simple_bandpass(temp, 2, 0);
    maxi=0;
    maxval=0;
    for i=packet_size+1:1000:length(temp)-packet_size-1
        if(sumsqr(real(temp(i:i+packet_size)))/sumsqr(real(temp(i-packet_size:i)))>maxval)
            maxi=i;
            maxval=sumsqr(real(temp(i:i+packet_size)))/sumsqr(real(temp(i-packet_size:i)));
        end
    end
    for i=maxi-1000:1:maxi-1000
        if(sumsqr(real(temp(i:i+packet_size)))/sumsqr(real(temp(i-packet_size:i)))>maxval)
            maxi=i;
            maxval=sumsqr(real(temp(i:i+packet_size)))/sumsqr(real(temp(i-packet_size:i)));
        end
    end
    rval1=maxi;
    temp=[temp(1:maxi-0.25*packet_size);temp(maxi+2*packet_size:end)];
    upspecial=simple_bandpass(up,2,r);
    downspecial=simple_bandpass(down,2,r);
    appSpecial=[downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;upspecial;upspecial;];
    maxi=zeros(round((length(temp)-length(downspecial))/1000),1);
    for i=1:1000:length(temp)-length(appSpecial)-1
        maxi(round(i/1000)+1)=max(abs(fftshift(fft(temp(i:i+length(appSpecial)-1).*appSpecial))));
    end
    if(rval2>rval1)
        [maxval maxloc]=max(maxi(1:end-(rval2-rval1)/1000));
    else
        [maxval maxloc]=max(maxi((rval1-rval2)/1000:end));
    end
    maxi=zeros(2001,1);
    maxloc=maxloc*1000+1;
    for i=maxloc-1000:maxloc+1000
        maxi(i-maxloc+1001)=max(abs(fftshift(fft(temp(i:i+length(appSpecial)-1).*appSpecial))));
    end
    [maxval1 maxloc1]=max(maxi);
    maxloc=maxloc-1001+maxloc1;
    eval((sprintf('s%d=temp;', 1)))
    ps_struct=[ps_struct;maxloc];
    loc_weight=[loc_weight;maxval1];
    
disp("checkpoint alpha")
toc
    
for rec=2:4
    data=rdata_struct{rec};
    temp=1i*data(:, 1)+data(:, 2);
    temp=simple_bandpass(temp, 2, r);
    temp=[temp(1:rval1-0.25*packet_size);temp(rval1+2*packet_size:end)];
    upspecial=simple_bandpass(up,2,r);
    downspecial=simple_bandpass(down,2,r);
    appSpecial=[downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;upspecial;upspecial;];
    maxi=zeros(round((length(temp)-length(downspecial))/1000),1);
    for i=1:1000:length(temp)-length(appSpecial)-1
        maxi(round(i/1000)+1)=max(abs(fftshift(fft(temp(i:i+length(appSpecial)-1).*appSpecial))));
    end
    [maxval maxloc]=max(maxi(round(ps_struct(1)/1000-packet_size/1000):round(ps_struct(1)/1000+packet_size/1000)));
    maxi=zeros(4001,1);
    maxloc=maxloc-1+round(ps_struct(1)/1000-packet_size/1000);
    maxloc=maxloc*1000+1;
    for i=maxloc-2000:maxloc+2000
        maxi(i-maxloc+2001)=max(abs(fftshift(fft(temp(i:i+length(appSpecial)-1).*appSpecial))));
    end
    [maxval1 maxloc1]=max(maxi);
    maxloc=maxloc-2001+maxloc1;
    eval((sprintf('s%d=temp;', rec)))
    ps_struct=[ps_struct;maxloc];
    loc_weight=[loc_weight;maxval1];
end

disp("checkpunt beta")
toc


data=rdata_struct{5};
    temp=1i*data(:, 1)+data(:, 2);
    temp=simple_bandpass(temp, 2, r);
    maxi=0;
    maxval=0;
    for i=packet_size+1:1000:length(temp)-packet_size-1
        if(sumsqr(real(temp(i:i+packet_size)))/sumsqr(real(temp(i-packet_size:i)))>maxval)
            maxi=i;
            maxval=sumsqr(real(temp(i:i+packet_size)))/sumsqr(real(temp(i-packet_size:i)));
        end
    end
    for i=maxi-1000:1:maxi-1000
        if(sumsqr(real(temp(i:i+packet_size)))/sumsqr(real(temp(i-packet_size:i)))>maxval)
            maxi=i;
            maxval=sumsqr(real(temp(i:i+packet_size)))/sumsqr(real(temp(i-packet_size:i)));
        end
    end
    rval2=maxi;
    temp=[temp(1:maxi-0.25*packet_size);temp(maxi+2*packet_size:end)];
    upspecial=simple_bandpass(up,2,r);
    downspecial=simple_bandpass(down,2,r);
    appSpecial=[downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;upspecial;upspecial;];
    maxi=zeros(round((length(temp)-length(downspecial))/1000),1);
    for i=1:1000:length(temp)-length(appSpecial)-1
        maxi(round(i/1000)+1)=max(abs(fftshift(fft(temp(i:i+length(appSpecial)-1).*appSpecial))));
    end
    [maxval maxloc]=max(maxi(round((ps_struct(1)+rval2-rval1)/1000-packet_size/1000):round((ps_struct(1)+rval2-rval1)/1000+packet_size/1000)));
    maxi=zeros(4001,1);
    maxloc=maxloc-1+round((ps_struct(1)+rval2-rval1)/1000-packet_size/1000);
    maxloc=maxloc*1000+1;
    for i=maxloc-2000:maxloc+2000
        maxi(i-maxloc+2001)=max(abs(fftshift(fft(temp(i:i+length(appSpecial)-1).*appSpecial))));
    end
    [maxval1 maxloc1]=max(maxi);
    maxloc=maxloc-2001+maxloc1;
    eval((sprintf('s%d=temp;', 5)))
    ps_struct=[ps_struct;maxloc];
    loc_weight=[loc_weight;maxval1];
    
disp("checkpoint charlie")
toc
    
for rec=6:8
    data=rdata_struct{rec};
    temp=1i*data(:, 1)+data(:, 2);
    temp=simple_bandpass(temp, 2, r);
    temp=[temp(1:rval2-0.25*packet_size);temp(rval2+2*packet_size:end)];
    upspecial=simple_bandpass(up,2,r);
    downspecial=simple_bandpass(down,2,r);
    appSpecial=[downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;upspecial;upspecial;];
    maxi=zeros(round((length(temp)-length(downspecial))/1000),1);
    for i=1:1000:length(temp)-length(appSpecial)-1
        maxi(round(i/1000)+1)=max(abs(fftshift(fft(temp(i:i+length(appSpecial)-1).*appSpecial))));
    end
    [maxval maxloc]=max(maxi(round(ps_struct(5)/1000-packet_size/1000):round(ps_struct(5)/1000+packet_size/1000)));
    maxi=zeros(4001,1);
    maxloc=maxloc-1+round(ps_struct(5)/1000-packet_size/1000);
    maxloc=maxloc*1000+1;
    for i=maxloc-2000:maxloc+2000
        maxi(i-maxloc+2001)=max(abs(fftshift(fft(temp(i:i+length(appSpecial)-1).*appSpecial))));
    end
    [maxval1 maxloc1]=max(maxi);
    maxloc=maxloc-2001+maxloc1;
    eval((sprintf('s%d=temp;', rec)))
    ps_struct=[ps_struct;maxloc];
    loc_weight=[loc_weight;maxval1];
end
    % rdata2=importdata('p1_sf7_2_IQ.csv');
% rdata3=importdata('13_3_IQ.csv');
% rdata4=importdata('13_4_IQ.csv');
% rdata5=importdata('13_5_IQ.csv');
% rdata6=importdata('13_6_IQ.csv');
% rdata7=importdata('13_7_IQ.csv');
% rdata8=importdata('13_8_IQ.csv');
 first_mean=mean(ps_struct(1:4))
  second_mean=mean(ps_struct(5:8))
  kark=round(second_mean-first_mean)-rval2+rval1;
  if kark<0
     kark=mod(kark,chirp_size)-chirp_size;
  else
      kark=mod(kark,chirp_size);
  end
  ps_struct(1:4)=repmat(round(first_mean+kark),4,1);
  ps_struct(5:8)=repmat(round(second_mean),4,1); 

  ps_struct(1:8)=ps_struct(1:8)+[0;ceil((length(s2)-length(s1))/2);ceil((length(s3)-length(s1))/2);ceil((length(s4)-length(s1))/2);0;ceil((length(s6)-length(s5))/2);ceil((length(s7)-length(s5))/2);ceil((length(s8)-length(s5))/2)];

disp("checkpoint delta")
toc
  
max_receiver=8;
p=zeros(chirp_size*packet_length,max_receiver);
n=zeros(chirp_size,max_receiver);
for i=1:max_receiver
    eval((sprintf('p(:,i)=s%d(ps_struct(i):ps_struct(i)+packet_size-1);', i)));
    eval((sprintf('n(:,i)=s%d(ps_struct(i)-5*chirp_size:ps_struct(i)-4*chirp_size-1);', i)));
end

   % %len=[length(rdata1) length(rdata2) length(rdata3) length(rdata4) length(rdata5) length(rdata6) length(rdata7) length(rdata8)];
% l=min(len);
% r=0;
% s=zeros(l,max_receiver);
% 
% 
% kk = 1;
% for counter=[1:8]
%     data=rdata_struct{counter};
%     jj=1;
%     nfactor=1;
%     nlocs=length(1:(chirp_size/nfactor):length(data)-nsym*chirp_size-1);
%     up_loc=zeros(nlocs, 1);
%     down_loc=zeros(nlocs, 1);
%     startloc=zeros(nlocs, 1);
%     for start_idx=1:(chirp_size/nfactor):length(data)-nsym*chirp_size-1
%         fourier_data=abs(fftshift(fft(data(start_idx:start_idx+chirp_size-1,1).*down)));
%         fourier_data_conj=abs(fftshift(fft(data(start_idx:start_idx+chirp_size-1,1).*up)));
%         [~, up_loc(jj)]=max(fourier_data);
%         [~, down_loc(jj)]=max(fourier_data_conj);
%         startloc(jj) = start_idx;
%         jj = jj+1;    
%     end
% 
%     likelihood = zeros(nlocs, 1);
%     for jj=1:nlocs-19*nfactor
%         likelihood(jj)= 1./(0.1+abs(std([up_loc(jj:nfactor:jj+15*nfactor);]))); % (down_loc(jj+16*nfactor)-0); (down_loc(jj+17*nfactor)-0); (down_loc([jj+18*nfactor jj+19*nfactor])-3)])));
%     end
% 
% 
%     x_axis=1:(chirp_size/nfactor):length(data)-nsym*chirp_size-1;
%     y_axis=likelihood*1000;
%     idx=find(y_axis==max(y_axis));
%     if(length(idx)>1)
%         idx=min(idx);
%     end
%     Xidx=x_axis(idx);
%     ps_struct(kk)=Xidx;
%     kk = kk + 1;
% end
% 
% noise_win1= ps_struct(1)-chirp_size:ps_struct(1)-1; % enter here
% noise_win2= ps_struct(5)-chirp_size:ps_struct(5)-1;      % enter here
% offset(1:4) =0;% floor((len(1:4) - len(1))./2);
% offset(5:8) =0;% floor((len(5:8) - len(5))./2);
%for rec=1:max_receiver
%    rdata=rdata_struct{rec}; 
%    if(rec<=4)
%        ps=ps_struct(1);
%        noise_win=noise_win1;
%        
%    else
%        ps=ps_struct(5);
%        noise_win=noise_win2;
%     
%    end
   
   
%end

reference=1;
%[~,reference]=max(loc_weight);
h=ones(1,max_receiver);
SNR_individual=zeros(1,max_receiver);
for count=1: max_receiver
    if (count~=reference)            
      angrat=angle(p(1:chirp_size,count)./p(1:chirp_size,reference));
      skip = 1;
      for skippow=1:13
          angrat=angle(p(3*chirp_size:12*chirp_size,count)./p(3*chirp_size:12*chirp_size,reference));
          slope = wrapToPi(angle(trimmean(exp(1i*(angrat(skip+1:end)-angrat(1:end-skip))), 0)));
          p(:,count)=p(:,count).*exp(-1j.*slope./skip.*(1:chirp_size*packet_length)');
          %n(:,count)=n(:,count).*exp(-1j.*slope./skip.*(-chirp_size+1:0)');
          
          skip = skip * 2;
      end
      %angrat=angle(p(1:chirpsize*10,count)./p(1:chirpsize*10,reference));
      %slope = wrapToPi(angle(trimmean(exp(1i*(angrat(1001:end)-angrat(1:end-1000))), 0)));
      %p(:,count)=p(:,count).*exp(-1j.*slope/1000.*(1:chirp_size*packet_length)');
      %n(:,count)=n(:,count).*exp(-1j.*slope/1000.*(-chirp_size:-1)');
      intercept = angle(mean(p(:,count)./p(:,reference)));
      h(1,count)=exp(-1j*intercept)*(mean(abs(p(:,count)./p(:,reference)).^2)^0.5);
      figure; hold on; title(count); plot(angle(h(:,count).*p(:,count)./p(:,reference)))
    end
    
    SNR_individual(1,count)=10*log10(max(abs(fftshift(fft(p(1:chirp_size,count).*down))))^2/mean(abs(fftshift(fft(n(1:chirp_size,count).*down))).^2));
end

disp("checkpoint echo")
toc

p_combined=p(:,1);
n_combined=n(:,1);
for ch=2:8
%       if (ch~=5 && ch~=6 && ch~=8)
    p_combined=p_combined+h(1,ch).*p(:,ch);
    n_combined=n_combined+n(:,ch);
%       end
end
% p_combined=sum(h.*n,2);
% n_combined=sum(h.*n,2);

SNR_combined=10*log10(max(abs(fftshift(fft(p_combined(1:chirp_size).*down))))^2/mean(abs(fftshift(fft(n_combined(1:chirp_size).*down))).^2));

disp("computed SNR combined")
toc