function [base_chirp,base_chirp_conj] = my_create_chirpspecialsmall(Fs,Ts,reset_freq,final_freq)
t=0:1/(Fs):Ts-(1/(Fs)); 

base_chirp=chirp(t,reset_freq/(2*2.315),t(end),final_freq/(2*2.315),'linear',90)+...
    1i*chirp(t,reset_freq/(2*2.315),t(end),final_freq/(2*2.315),'linear'); 
base_chirp=reshape(base_chirp,length(base_chirp),1);
base_chirp_conj = conj(base_chirp);
end