function [base_chirp,base_chirp_conj] = my_create_chirpspecial1(Fs,Ts,reset_freq,final_freq,symbol_length)
    t=0:1/(Fs):Ts-(1/(Fs)); 

    base_chirp = chirp(t, reset_freq*1024/symbol_length, t(end),final_freq * 1024 / symbol_length, 'linear', 90) + 1i*chirp(t,reset_freq*1024/symbol_length,t(end),final_freq*1024/symbol_length,'linear'); 

    base_chirp = reshape(base_chirp,length(base_chirp),1);
    base_chirp_conj = conj(base_chirp);
end