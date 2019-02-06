function [base_chirp] = new_chirp(Fs,Ts,reset_freq,final_freq,symbol_length)

Tsamp = 1/Fs; % sampling frequency?
t=0:Tsamp:Ts-Tsamp; 
tprime=round(rand()*length(t));
slope=(final_freq-reset_freq)*1024/symbol_length;
ender=(t(end)-tprime)*slope+reset_freq*1024/symbol_length;
base_chirp1=[chirp(t,ender,tprime,final_freq*1024/symbol_length,'linear',90)]+...
    1j*[chirp(t,ender,tprime,final_freq*1024/symbol_length,'linear',90) ];
base_chirp2= chirp(tprime,reset_freq*1024/symbol_length,t(end),ender,'linear',90) +1j*chirp(tprime,reset_freq*1024/symbol_length,t(end),ender,'linear',90);
base_chirp=complex(zeros(length(t),1));
base_chirp(1:tprime)=base_chirp(1:tprime)+transpose(base_chirp1);
base_chirp(tprime:end)=base_chirp(tprime:end)+transpose(base_chirp2);

end