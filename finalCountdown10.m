%% setup

BW = 1.25e5;
SF = 10;
r = 0;
chirp_size = 4608; % length of the chirp in number of samples
packet_length = 20; % number of chirps in a packet
packet_size = chirp_size*packet_length;
bandwidth_sampling_factor = 1;

% table  Sample/KHz       bandwidth_sampling_factor
%         125                1
%         250                2
%         500                3

Fs = BW;
symbol_length = chirp_size;
symbol_length_upsampled = bandwidth_sampling_factor*chirp_size;
freq_shift_per_sample =  Fs/symbol_length; % How each frequency bin maps to a difference in frequency
Ts = 1/freq_shift_per_sample; % Symbol Duration
f = linspace(-BW/2,BW/2-freq_shift_per_sample,symbol_length); % The X-Axis of the FFT plot
reset_freq = -BW/2; % The initial frequency of the base chirp
final_freq = (BW/2);%-freq_shift_per_sample; % The final frequency
nsym = 38.25;
len = zeros(1, 8);
rdata_struct = cell(1, 8);

%% data import

% import data from various gateways
for ii = 1:8
	rdata_struct{ii} = importdata(sprintf('data/p2_sf10_%d_IQ.csv', ii));
	len(ii) = length(rdata_struct{ii});
end

%%% why is detection being done seperately?
% packet detection on stream 5

data = rdata_struct{5};
temp = 1i*data(:, 1) + data(:, 2); % convert IQ samples to complex numbers
temp = simple_bandpass(temp, 2, r); % bandpass filter???

maxi = 0;
maxval = 0;

% windowed energy detection to find the packet
for i = (packet_size + 1):1000:(length(temp) - packet_size - 1)
	if(sumsqr(real(temp(i:i+packet_size)))/sumsqr(real(temp(i-packet_size:i)))>maxval)
		maxi = i;
		maxval = sumsqr(real(temp(i:i+packet_size)))/sumsqr(real(temp(i-packet_size:i)));
	end
end
rval2 = maxi;

% packet detection on stream 1

data = rdata_struct{1};
temp = 1i*data(:, 1) + data(:, 2);
temp = simple_bandpass(temp, 2, r);

maxi = 0;
maxval = 0;

for i = (packet_size + 1):1000:(length(temp) - packet_size - 1)
	if(sumsqr(real(temp(i:i+packet_size)))/sumsqr(real(temp(i-packet_size:i)))>maxval)
		maxi = i;
		maxval = sumsqr(real(temp(i:i+packet_size)))/sumsqr(real(temp(i-packet_size:i)));
	end
end

rval1 = maxi;


temp = [temp(1:(maxi - 0.25 * packet_size)); temp((maxi + 2 * packet_size):end)];

%%
% create base upchirp and downchirps

[up,down] = my_create_chirpspecial1(bandwidth_sampling_factor*Fs,Ts,reset_freq,final_freq,chirp_size);

upspecial = simple_bandpass(up,2,r);
downspecial = simple_bandpass(down,2,r);

% WHAT IS THIS?
% 18 downchirps and 2 upchirps
appSpecial = [downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;upspecial;upspecial;];

maxi = zeros(round((length(temp)-length(downspecial))/1000),1);
for i = 1:1000:length(temp)-length(appSpecial)-1
	maxi(round(i/1000)+1) = max(abs(fftshift(fft(temp(i:i+length(appSpecial)-1).*appSpecial))));
end

if(rval2>rval1)
	[maxval maxloc]=max(maxi(1:end-(rval2-rval1)/1000));
	maxloc=maxloc*1000+1;
else
	[maxval maxloc]=max(maxi((rval1-rval2)/1000:end));
	maxloc=(maxloc+(rval1-rval2)/1000)*1000+1;

end
maxi=zeros(2001,1);
for i=maxloc-1000:maxloc+1000
	maxi(i-maxloc+1001)=max(abs(fftshift(fft(temp(i:i+length(appSpecial)-1).*appSpecial))));
end
[maxval1 maxloc1]=max(maxi);
maxloc=maxloc-1001+maxloc1;
eval((sprintf('s%d=temp;', 1)))

ps_struct = [];
loc_weight = [];

ps_struct=[ps_struct;maxloc];
loc_weight=[loc_weight;maxval1];

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
	[maxval maxloc]=max(maxi(round(ps_struct(1)/1000-packet_size/8000):round(ps_struct(1)/1000+packet_size/8000)));
	maxi=zeros(4001,1);
	maxloc=maxloc-1+round(ps_struct(1)/1000-packet_size/8000);
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
  %ps_struct(1:4)=repmat(round(first_mean),4,1);
  ps_struct(5:8)=repmat(round(second_mean),4,1); 

  ps_struct(1:8)=ps_struct(1:8)+[0;floor((length(s2)-length(s1))/2);floor((length(s3)-length(s1))/2);floor((length(s4)-length(s1))/2);0;floor((length(s6)-length(s5))/2);floor((length(s7)-length(s5))/2);floor((length(s8)-length(s5))/2)];


max_receiver=8;
p=zeros(chirp_size*packet_length,max_receiver);
n=zeros(chirp_size,max_receiver);
for i=1:max_receiver
	eval((sprintf('p(:,i)=s%d(ps_struct(i):ps_struct(i)+packet_size-1);', i)));
	eval((sprintf('n(:,i)=s%d(ps_struct(i)-chirp_size:ps_struct(i)-1);', i)));
end
   
   
%end

reference=1;
%[~,reference]=max(loc_weight);
h=ones(1,max_receiver);
SNR_individual=zeros(1,max_receiver);

for count=1: max_receiver
	if (count~=reference)            
	  angrat=angle(p(1:chirp_size,count)./p(1:chirp_size,reference));
	  skip = 1;
	  for skippow=1:16
		  angrat=angle(p(3*chirp_size:12*chirp_size,count)./p(3*chirp_size:12*chirp_size,reference));
		  slope = wrapToPi(angle(mean(exp(1j*(angrat(skip+1:end)-angrat(1:end-skip))))));
		  p(:,count)=p(:,count).*exp(-1j.*slope./skip.*(1:chirp_size*packet_length)');
		  n(:,count)=n(:,count).*exp(-1j.*slope./skip.*(-chirp_size+1:0)');
		  skip = skip * 2;
	  end
	  %angrat=angle(p(1:chirpsize*10,count)./p(1:chirpsize*10,reference));
	  %slope = wrapToPi(angle(trimmean(exp(1i*(angrat(1001:end)-angrat(1:end-1000))), 0)));
	  %p(:,count)=p(:,count).*exp(-1j.*slope/1000.*(1:chirp_size*packet_length)');
	  %n(:,count)=n(:,count).*exp(-1j.*slope/1000.*(-chirp_size:-1)');
	  intercept = angle(mean(p(:,count)./p(:,reference)));
	  h(1,count) = exp(-1j*intercept)*(mean(abs(p(:,count)./p(:,reference)).^2)^0.5);
	  figure;
	  hold on;
	  title(count);
	  plot(angle(h(:,count).*p(:,count)./p(:,reference)))
	end
	
	SNR_individual(1,count)=10*log10(max(abs(fftshift(fft(p(1:chirp_size,count).*down))))^2/mean(abs(fftshift(fft(n(1:chirp_size,count).*down))).^2));
end

% compute combined signal

p_combined=p(:,1);
n_combined=n(:,1);
for ch=2:8
%       if (ch~=5 && ch~=6 && ch~=8)
	p_combined=p_combined+h(1,ch).*p(:,ch);
	n_combined=n_combined+h(1,ch).*n(:,ch);
%       end
end
% p_combined=sum(h.*n,2);
% n_combined=sum(h.*n,2);

SNR_combined=10*log10(max(abs(fftshift(fft(p_combined(1:chirp_size).*down))))^2/mean(abs(fftshift(fft(n_combined(1:chirp_size).*down))).^2));