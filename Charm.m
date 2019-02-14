SF = 10;
BW = 1.25e5;
chirp_size=1024;
chirpsize=chirp_size;
extra_sampling_factor = 1;
Fs = BW;
symbol_length=chirp_size;
symbol_length_upsampled = extra_sampling_factor*chirp_size;
freq_shift_per_sample =  Fs/symbol_length; % How each frequency bin maps to a difference in frequency
Ts = 1/freq_shift_per_sample; % Symbol Duration
f = linspace(-BW/2,BW/2-freq_shift_per_sample,symbol_length); % The X-Axis
reset_freq = -BW/2; % The initial frequency of the base chirp
final_freq = (BW/2)-freq_shift_per_sample; % The final frequency
[up,down] =  my_create_chirpspecial(extra_sampling_factor*Fs,Ts,reset_freq,final_freq,chirp_size);
upfft=fft(up);
downfft=fft(down);
up250=(ifft([upfft(1:512);zeros(1024,1);upfft(513:1024)]));
down250=(ifft([downfft(1:512);zeros(1024,1);downfft(513:1024)]));
sig1=read_complex_binary('predabs52_15_04_47.dat');
sig2=read_complex_binary('predabs3_15_04_47.dat');
sig3=read_complex_binary('predabs7_15_04_47.dat');
sig4=read_complex_binary('predabs50_15_04_47.dat');

p1=getPacket(sig1,down250,up250);
p2=getPacket(sig2,down250,up250);
p3=getPacket(sig3,down250,up250);
p4=getPacket(sig4,down250,up250);
close all;
figure; hold on;
plot(abs(p1+p2+p3+p4),'Color','black');
plot(abs(p1),'Color','blue')
plot(abs(p2),'Color','blue')
plot(abs(p3),'Color','blue')
plot(abs(p4),'Color','blue')

p2=OffsetCorrectorNew(p2,p1);
p3=OffsetCorrectorNew(p3,p1);
p4=OffsetCorrectorNew(p4,p1);

plot(abs(p1+p2+p3+p4),'Color','red')


