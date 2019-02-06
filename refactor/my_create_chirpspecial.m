function [base_chirp,base_chirp_conj] = my_create_chirpspecial(Fs,Ts,reset_freq,final_freq)
Tsample = 1/Fs; % sampling frequency
t=0:Tsample:Ts-Tsample; 

base_chirp=chirp(t,reset_freq/(2*2.315),t(end),final_freq/(2*2.315),'linear',90)+...
    1i*chirp(t,reset_freq/(2*2.315),t(end),final_freq/(2*2.315),'linear'); 
base_chirp=reshape(base_chirp,length(base_chirp),1);
base_chirp_conj = conj(base_chirp);
end