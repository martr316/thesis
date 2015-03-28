function [fit, RMSE  ] = validate_model(model, measured_data, Fs)
%VALIDATE_MODEL simulates the model with measured data and plots the output
%in amplitude and frequency. The output from function is the
%fit-value from the matlabfunction 'compare' and a RMSE value from the error in
%amplitude.

load(measured_data);
loaded_model=load(model);
fields=fieldnames(loaded_model);

Total_samples = MeasureddataExcersion.Total_Samples;
time=(1:Total_samples)/Fs;

% simulate and calculate resonse and fit
[simulated_data,~]=lsim(loaded_model.(fields{1}),MeasureddataVoltage.Data,time);
data = iddata(MeasureddataExcersion.Data/1000,MeasureddataVoltage.Data,1/Fs); %divide by 1000 to get in SI-units
[~,fit,~] = compare(data,loaded_model.(fields{1})); %calculate the response and the fit 

figure(1)
plot(time,simulated_data,time,MeasureddataExcersion.Data/1000)
legend('Simulated data','Measured data')

%calculate root mean square error for amplitude
RMSE = sqrt(mean(MeasureddataExcersion.Data/1000-simulated_data).^2);

% Values for fft-calcualtions
NFFT = 2^nextpow2(Total_samples); % Next power of 2 from length of y
f = Fs/2*linspace(0,1,NFFT/2+1);

%Find the frequency resonse from model
simulated_fft = fft(simulated_data(1:Total_samples),NFFT)/Total_samples;
db_simulated_fft = mag2db(2*abs(simulated_fft(1:NFFT/2+1))); %convert to decibel scale

%Find the frequency response from measured data
measured_fft = fft(MeasureddataExcersion.Data(1:Total_samples),NFFT)/Total_samples;
db_measured_fft = mag2db(2*abs(measured_fft(1:NFFT/2+1))); %convert to decibel scale

% PLot FFT
figure(2)
loglog(f,db_simulated_fft,f,db_measured_fft)
title('FFT of validation data')
xlabel('Frequency (Hz)')
ylabel('|Input(f)|')
legend('Simulated data','Measured data')

% %Find the frequency resonse from model
% [mag,ph,w]=bode(H_total,20*(2*pi):20000*(2*pi)); %from 20Hz-20kHz
% mag=squeeze(mag);
% plot(w/(2*pi),mag2db(mag))
% hold off
end

