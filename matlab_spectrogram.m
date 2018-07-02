filename = "data/914.5MHz_2048kHz_8dB_cap";
file = fopen(filename, 'r');
data = fread(file, [2, Inf], 'float32').';
fclose(file);

data = data(:,1) + 1i * data(:,2); 
disp([num2str(size(data, 1)), ' samples read'])
fs = 2048000;

%%
figure(1);
subplot(2,1,1);
window_size = 512;
spectrogram(data(1:204800*1), hann(window_size), window_size./2, window_size, 'yaxis');

subplot(2,1,2);
spectrogram(data(204800*10:204800*11), hann(window_size), window_size./2, window_size, 'yaxis');

