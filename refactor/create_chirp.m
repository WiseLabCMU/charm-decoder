function [base_chirp] = create_chirp(Fs,Ts,reset_freq,final_freq,symbol_length)
    
    t = 0:1/(Fs):Ts-(1/(Fs)); 
    Fstart = reset_freq*1024/symbol_length;
    Fstop = final_freq*1024/symbol_length;
    I = chirp(t, Fstart, t(end), Fstop, 'linear', 90);
    Q = chirp(t,Fstart,t(end),Fstop,'linear');
    base_chirp = I + 1i*Q;
    base_chirp = reshape(base_chirp,length(base_chirp),1); % return column vector
end