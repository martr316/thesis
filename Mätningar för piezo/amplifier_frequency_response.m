
%%
%a = simpleConvertTDMS('förstärkare_kapacitans_40kHz_29nF.tdms');
load('förstärkare_kapacitans_40kHz_29nF')
data2=MeasureddataVoltage.Data;

Fs=44100;
L = length(MeasureddataVoltage.Data);

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
out_fft = fft(MeasureddataVoltage.Data,NFFT)/L;
in_fft = fft(MeasureddataInput.Data,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

amplifier_fft = out_fft./in_fft;
fft_freq =2*abs(amplifier_fft(1:NFFT/2+1));
% Plot single-sided amplitude spectrum.
figure(57)
plot(f,mag2db(fft_freq))
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
hold on


cutoff_sample=12e4;
p = polyfit(1:cutoff_sample,fft_freq(1:cutoff_sample)',2);
y1=polyval(p,1:cutoff_sample);
plot(f(1:cutoff_sample),mag2db(y1))
hold off

%%
load('förstärkare_kapacitans_29nF')
Total_samples = MeasureddataExcersion.Total_Samples;
data1=MeasureddataVoltage.Data;
figure(1)
plot(abs(fft(data1)));

out=fft(MeasureddataVoltage.Data);
in =fft(MeasureddataInput.Data);
BB= out./in;

figure(58)
plot(mag2db(abs(BB)))

%%
a = simpleConvertTDMS('förstärkare_resistans_40kHz_68ohm.tdms');
load('förstärkare_resistans_40kHz_68ohm')
data3=MeasureddataVoltage.Data;
figure(3)
plot(abs(fft(data2)));


out=fft(MeasureddataVoltage.Data);
in =fft(MeasureddataInput.Data);
BB= out./in;

figure(57)
plot(mag2db(abs(BB)))
