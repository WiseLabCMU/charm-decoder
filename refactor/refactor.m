tic
% Transmitted chirp and packet parameters
BW = 125000.0; % chirp bandwidth is 125 kHz
SF = 7; % LoRa spreading factor of the transmission
packet_length = 20; % size of the packet in number of chirps (?)
symbolTime = pow2(SF)/BW; % from LoRa manual


% receiver parameters
f_xo = 36e6; % frequency of receiver's crystal oscillator
sigmaDeltaWindow = 64; % number of 1 bit samples to combine for a sample value
decimationFactor = 4; % reducing our sampling from 562.5 kHz to 140.625 kHz
% this was set in the LPRAN FPGA code.

samplesPerSec = f_xo/(sigmaDeltaWindow);


% potential problem: we seem to use the samplesPerSymbol value assuming
% ~500 kHz bandwidth but use other parameters (upsampling factor) from assuming ~125 kHz
samplesPerSymbol = round(symbolTime * samplesPerSec); % This should be 576 samples for 125 kHz
packet_size = samplesPerSymbol*packet_length;

% table  Sample/KHz       extra_sampling_factor
%         125                1
%         250                2
%         500                3
% this functions should just be factor = 1 + log2(samplingFreq/bandwidth);
% TODO: why is the factor for 500 kHz = 3 and not 4?

% TODO: here is a big issue to figure out. The reference chirps seem to be
% generated assuming a 125 kHz bandwidth while the actual signal seems to
% be captured with 500 kHz of bandwidth;

upsamplingFactor = round(1 + log2(samplesPerSec/(decimationFactor*BW)));
% upsamplingFactor = 1;
upsampledSamplesPerSymbol = upsamplingFactor*samplesPerSymbol;

freq_shift_per_sample =  BW/samplesPerSymbol; % How each frequency bin maps to a difference in frequency

% f = linspace(-BW/2, BW/2-freq_shift_per_sample, samplesPerSymbol); % The X-Axis of the FFT plot

nsym = 38.25;

disp("initialized parameters")
toc

% import data from .csv files
% these need to be cells as the different files may have different lengths
len = zeros(1, 8);
rdata = cell(1, 8);
for index=1:8
    temp = importdata(sprintf('data/p1_sf7_%d_IQ.csv', index));
    rdata{index} = 1i*temp(:,1) + temp(:,2);
end

disp("loaded data")
toc

% process the 5th reception (picked at random)
randomSampleIndex = 5; % pick a random reception as our reference
temp = simple_bandpass(rdata{randomSampleIndex}, 2, 0);
maxi = 0;
maxval = 0;
% energy detection to find the preamble
% Note: this would get screwy if there is more than one packet in the window
% Note 2: This could be easily parallelized


% first do a coarse search for every 1000th sample
for i = packet_size+1:1000:length(temp)-packet_size-1
    
    % Q: Why are only the real signals (I channel) looked at for energy
    % detection
    energyNext = sumsqr(real(temp(i:i+packet_size)));
    energyPrev = sumsqr(real(temp(i-packet_size:i)));
    if(energyNext/energyPrev > maxval)
        maxi=i;
        maxval = energyNext/energyPrev;
    end
end

% then do a finer grained search for the largest candidate
for i=maxi-1000:1:maxi-1000 % this seems to be wrong i should vary from maxi-1000 to maxi
    energyNext = sumsqr(real(temp(i:i+packet_size)));
    energyPrev = sumsqr(real(temp(i-packet_size:i)));
    if(energyNext/energyPrev > maxval)
        maxi=i;
        maxval = energyNext/energyPrev;
    end
end
rval2=maxi; % record the location of the peak

loc_weight=[];
randomSampleIndex = 1; % pick a random reception as our reference
temp = simple_bandpass(rdata{randomSampleIndex}, 2, 0);
maxi=0;
maxval=0;

% coarse grained search for another random sample
for i=packet_size+1:1000:length(temp)-packet_size-1
    energyNext = sumsqr(real(temp(i:i+packet_size)));
    energyPrev = sumsqr(real(temp(i-packet_size:i)));
    if(energyNext/energyPrev > maxval)
        maxi=i;
        maxval = energyNext/energyPrev;
    end
end

% fine grained search for second random sample
for i=maxi-1000:1:maxi-1000 % this seems to be wrong i should vary from maxi-1000 to maxi
    energyNext = sumsqr(real(temp(i:i+packet_size)));
    energyPrev = sumsqr(real(temp(i-packet_size:i)));
    if(energyNext/energyPrev > maxval)
        maxi=i;
        maxval = energyNext/energyPrev;
    end
end
rval1 = maxi; % record the location of the peak

% extract a smaller chunk out of the buffer
% this one only seems to consider the region which definitely does not
% contain the packet
% NOTE: This might be for noise level calculation
disp(["samples taken", 1, "to", rval1-0.25*packet_size, "and", rval1+2*packet_size, "to", length(temp)])
temp = [temp(1:rval1-0.25*packet_size);temp(rval1+2*packet_size:end)];

disp("analyzed both random samples")
toc

% create reference chirps
reset_freq = -BW/2; % The initial frequency of the base chirp
final_freq = (BW/2);%-freq_shift_per_sample; % The final frequency
Ts = 1/freq_shift_per_sample;

up = create_chirp(upsamplingFactor*BW,Ts,reset_freq,final_freq,samplesPerSymbol); % perfect upchirp
down = conj(up); % perfect downchirp
upspecial = simple_bandpass(up,2,0); % band limited upchirp
downspecial = simple_bandpass(down,2,0); % band limited downchirp

% construct conjugate of expected pramble
appSpecial = [downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;upspecial;upspecial;];

% Q. is this the noise level calculation?
maxiSize = round((length(temp)-length(downspecial))/1000);
maxi = zeros(maxiSize,1);

% start with a coarse grained calculation
for i = 1:1000:length(temp)-length(appSpecial)-1
    filtered = temp(i:i+length(appSpecial)-1).*appSpecial;
    maxi(round(i/1000)+1) = max(abs(fftshift(fft( filtered ))));
end

% huh?
if(rval2 > rval1)
    [~, maxloc] = max(maxi(1:end-(rval2-rval1)/1000));
else
    [~, maxloc] = max(maxi((rval1-rval2)/1000:end));
end
maxi = zeros(2001,1);
maxloc = maxloc*1000+1;

% then a fine grained calculation
for i = maxloc-1000:maxloc+1000
    filtered = temp(i:i+length(appSpecial)-1).*appSpecial;
    maxi(i-maxloc+1001) = max(abs(fftshift(fft( filtered ))));
end
[maxval1, maxloc1] = max(maxi);
maxloc = maxloc - 1001 + maxloc1;


s = cell(1,8);
s{1} = temp;
% eval((sprintf('s%d=temp;', 1)))

disp("generated reference chirps")
toc

ps_struct = []; % what is this for???
ps_struct = [ps_struct;maxloc];
loc_weight=[loc_weight;maxval1];


for rec = 2:4
    temp = simple_bandpass(rdata{rec}, 2, 0);
    temp = [temp(1:rval1-0.25*packet_size); temp(rval1+2*packet_size:end)];
    
    
    upspecial=simple_bandpass(up,2,0);
    downspecial=simple_bandpass(down,2,0);
    
    % construct conjugate of preamble (for multiplication with original signal)
    appSpecial=[downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;upspecial;upspecial;];
    maxi=zeros(round((length(temp)-length(downspecial))/1000),1);
    for i=1:1000:length(temp)-length(appSpecial)-1
        maxi(round(i/1000)+1)=max(abs(fftshift(fft(temp(i:i+length(appSpecial)-1).*appSpecial))));
    end
    [~, maxloc]=max(maxi(round(ps_struct(1)/1000-packet_size/1000):round(ps_struct(1)/1000+packet_size/1000)));
    maxi=zeros(4001,1);
    maxloc=maxloc-1+round(ps_struct(1)/1000-packet_size/1000);
    maxloc=maxloc*1000+1;
    for i=maxloc-2000:maxloc+2000
        maxi(i-maxloc+2001)=max(abs(fftshift(fft(temp(i:i+length(appSpecial)-1).*appSpecial))));
    end
    [maxval1, maxloc1]=max(maxi);
    maxloc=maxloc-2001+maxloc1;
    s{rec} = temp;
    ps_struct=[ps_struct;maxloc];
    loc_weight=[loc_weight;maxval1];
end

disp("processed receptions 2 to 4")
toc

randomSampleIndex = 5; % pick a random reception as our reference
temp=simple_bandpass(rdata{randomSampleIndex}, 2, 0);
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
upspecial=simple_bandpass(up,2,0);
downspecial=simple_bandpass(down,2,0);
appSpecial=[downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;upspecial;upspecial;];
maxi=zeros(round((length(temp)-length(downspecial))/1000),1);
for i=1:1000:length(temp)-length(appSpecial)-1
    maxi(round(i/1000)+1)=max(abs(fftshift(fft(temp(i:i+length(appSpecial)-1).*appSpecial))));
end
[maxval, maxloc]=max(maxi(round((ps_struct(1)+rval2-rval1)/1000-packet_size/1000):round((ps_struct(1)+rval2-rval1)/1000+packet_size/1000)));
maxi=zeros(4001,1);
maxloc=maxloc-1+round((ps_struct(1)+rval2-rval1)/1000-packet_size/1000);
maxloc=maxloc*1000+1;
for i=maxloc-2000:maxloc+2000
    maxi(i-maxloc+2001)=max(abs(fftshift(fft(temp(i:i+length(appSpecial)-1).*appSpecial))));
end
[maxval1, maxloc1]=max(maxi);
maxloc=maxloc-2001+maxloc1;
% eval((sprintf('s%d=temp;', 5)))
s{5} = temp;
ps_struct=[ps_struct;maxloc];
loc_weight=[loc_weight;maxval1];

for rec=6:8
    temp=simple_bandpass(rdata{rec}, 2, 0);
    temp=[temp(1:rval2-0.25*packet_size);temp(rval2+2*packet_size:end)];
    upspecial=simple_bandpass(up,2,0);
    downspecial=simple_bandpass(down,2,0);
    appSpecial=[downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;downspecial;upspecial;upspecial;];
    maxi=zeros(round((length(temp)-length(downspecial))/1000),1);
    for i=1:1000:length(temp)-length(appSpecial)-1
        maxi(round(i/1000)+1)=max(abs(fftshift(fft(temp(i:i+length(appSpecial)-1).*appSpecial))));
    end
    [maxval, maxloc]=max(maxi(round(ps_struct(5)/1000-packet_size/1000):round(ps_struct(5)/1000+packet_size/1000)));
    maxi=zeros(4001,1);
    maxloc=maxloc-1+round(ps_struct(5)/1000-packet_size/1000);
    maxloc=maxloc*1000+1;
    for i=maxloc-2000:maxloc+2000
        maxi(i-maxloc+2001)=max(abs(fftshift(fft(temp(i:i+length(appSpecial)-1).*appSpecial))));
    end
    [maxval1, maxloc1]=max(maxi);
    maxloc=maxloc-2001+maxloc1;
%     eval((sprintf('s%d=temp;', rec)))
    s{rec} = temp;
    ps_struct=[ps_struct;maxloc];
    loc_weight=[loc_weight;maxval1];
end

disp("processed the other receptions")
toc

first_mean = mean(ps_struct(1:4));
second_mean = mean(ps_struct(5:8));
kark=round(second_mean-first_mean)-rval2+rval1;
if kark<0
    kark=mod(kark,samplesPerSymbol)-samplesPerSymbol;
else
    kark=mod(kark,samplesPerSymbol);
end
ps_struct(1:4)=repmat(round(first_mean+kark),4,1);
ps_struct(5:8)=repmat(round(second_mean),4,1);

ps_struct(1:8)=ps_struct(1:8) + ...
        [   
            0; ...
            ceil((length(s{2})-length(s{1}))/2); ...
            ceil((length(s{3})-length(s{1}))/2); ...
            ceil((length(s{4})-length(s{1}))/2); ...
            0; ...
            ceil((length(s{6})-length(s{5}))/2); ...
            ceil((length(s{7})-length(s{5}))/2); ...
            ceil((length(s{8})-length(s{5}))/2)
        ];


max_receiver=8;
p=zeros(samplesPerSymbol*packet_length,max_receiver);
n=zeros(samplesPerSymbol,max_receiver);
for i=1:max_receiver
    p(:,i)=s{i}(ps_struct(i):ps_struct(i)+packet_size-1);
    n(:,i)=s{i}(ps_struct(i)-5*samplesPerSymbol:ps_struct(i)-4*samplesPerSymbol-1);
end

reference=1;
h = ones(1,max_receiver);
SNR_individual=zeros(1,max_receiver);
for count=1: max_receiver
    if (count~=reference)
        angrat=angle(p(1:samplesPerSymbol,count)./p(1:samplesPerSymbol,reference));
        skip = 1;
        for skippow=1:13
            angrat=angle(p(3*samplesPerSymbol:12*samplesPerSymbol,count)./p(3*samplesPerSymbol:12*samplesPerSymbol,reference));
            slope = wrapToPi(angle(trimmean(exp(1i*(angrat(skip+1:end)-angrat(1:end-skip))), 0)));
            p(:,count)=p(:,count).*exp(-1j.*slope./skip.*(1:samplesPerSymbol*packet_length)');
            skip = skip * 2;
        end

        intercept = angle(mean(p(:,count)./p(:,reference)));
        h(1,count)=exp(-1j*intercept)*(mean(abs(p(:,count)./p(:,reference)).^2)^0.5);
        figure;
        hold on;
        title(count);
        plot(angle(h(:,count).*p(:,count)./p(:,reference)));
    end
    x = max(abs(fftshift(fft(p(1:samplesPerSymbol,count).*down))));
    size(x);
    y = mean(abs(fftshift(fft(n(1:samplesPerSymbol,count).*down))).^2);
    size(y);
    SNR_individual(1,count) = 10*log10(x^2/y);
end

p_combined = p(:,1);
n_combined = n(:,1);
for ch=2:8
    p_combined = p_combined+h(1,ch).*p(:,ch);
    n_combined = n_combined+n(:,ch);
end

x = max(abs(fftshift(fft(p_combined(1:samplesPerSymbol).*down))));
y = mean(abs(fftshift(fft(n_combined(1:samplesPerSymbol).*down))).^2);
SNR_combined=10*log10(x^2/y);