%% Find the damping constants for different modes
m = 2.6657e-04;
l = 0.0267;
Fs=44100;

%% First mode at 205hz 
a = simpleConvertTDMS('sinus_mode1_205hz.tdms');
load('sinus_mode1_205hz');

Total_samples = MeasureddataExcersion.Total_Samples;

plot((1:Total_samples)/Fs,MeasureddataExcersion.Data)
title('Excursion')
xlabel('Seconds')
ylabel('Milimeter')

fi1 = 0.02549;
fi2 = 0.01466;
tau1 = 3.012; %[s]
tau2 = 3.083; %[s]

q=15; %periods
tau_d_m= tau2-tau1; %[s]

r1 = (2 * m)/tau_d_m * log(fi1/fi2); %[Ns/m]
r_a1 = r1/l;

%% Second mode at 1430hz 
a = simpleConvertTDMS('sinus_mode2_1430hz.tdms');
load('sinus_mode2_1430hz');

Total_samples = MeasureddataExcersion.Total_Samples;

plot((1:Total_samples)/Fs,MeasureddataExcersion.Data)
title('Excursion')
xlabel('Seconds')
ylabel('Milimeter')

% fi1 = 0.005145;
% fi2 = 0.00007926;
% tau1 = 3.012; %[s]
% tau2 = 3.05; %[s]

fi1 = 0.02135;
fi2 = 0.004891;
tau1 = 3.007; %[s]
tau2 = 3.009; %[s]

q=3; %periods
tau_d_m= tau2-tau1; %[s]

r2 = (2 * m)/tau_d_m * log(fi1/fi2); %[Ns/m]
r_a2 = r2/l;

%% Third mode at 4300hz 
a = simpleConvertTDMS('sinus_mode3_4300hz.tdms');
load('sinus_mode3_4300hz');

Total_samples = MeasureddataExcersion.Total_Samples;

plot((1:Total_samples)/Fs,MeasureddataExcersion.Data)
title('Excursion')
xlabel('Seconds')
ylabel('Milimeter')

fi1 = 0.002036-0.0006551;
fi2 = 0.000972-0.000373;
tau1 = 3.008; %[s]
tau2 = 3.01; %[s]

q=2; %periods
tau_d_m= tau2-tau1; %[s]

r3 = (2 * m)/tau_d_m * log(fi1/fi2); %[Ns/m]
r_a3 = r3/l;


r = [ r1 r2 r3]
r_a = [ r_a1 r_a2 r_a3]


%% DO NOT USE!!! old test values
% Find the r_a constant for piezo

a = simpleConvertTDMS('recording1000Hz.tdms');

load('recording1000Hz.mat')
Fs=44100; % Sample frequency (Hz)

Total_samples = MeasureddataExcersion.Total_Samples;

plot(((Fs*0.9):(1.2*Fs))/Fs,MeasureddataExcersion.Data((Fs*0.9):(1.2*Fs)))
title('Excursion')
xlabel('Seconds')
ylabel('Milimeter')

% Values from plot
fi1 = 0.02099e-3;
fi2 = 0.001417e-3;
tau1 = 1.015; %[s]      %4.474e4; Samples
tau2 = 1.099; %[s]      %4.845e4;

q=18; %periods
tau_d_m= tau2-tau1; %[s]

m = 2.6657e-04;
l = 0.0267;

r = (2 * m)/tau_d_m * log(fi1/fi2) %[Ns/m]

r_a = r/l %[Ns/m^2] 

%% Test2 1000Hz
clear all
a = simpleConvertTDMS('recording1000Hz.tdms');

load('recording1000Hz.mat')
Fs=44100; % Sample frequency (Hz)

Total_samples = MeasureddataExcersion.Total_Samples;

plot((1:Total_samples)/Fs,MeasureddataExcersion.Data)
title('Excursion 1000Hz')
xlabel('Seconds')
ylabel('Milimeter')

% Values from plot
fi1 = 0.02114e-3;
fi2 = 0.001072e-3;
tau1 = 1.014; %[s]      %4.474e4; Samples
tau2 = 1.108; %[s]

q=20; %periods
tau_d_m= tau2-tau1; %[s]
m = 2.6657e-04;
l = 0.0267;
r = (2 * m)/tau_d_m * log(fi1/fi2) %[Ns/m]

r_a = r/l %[Ns/m^2] 

%% Test3 100Hz
clear all
a = simpleConvertTDMS('recording100Hz.tdms');

load('recording100Hz.mat')
Fs=44100; % Sample frequency (Hz)

Total_samples = MeasureddataExcersion.Total_Samples;

plot((1:Total_samples)/Fs,MeasureddataExcersion.Data)
title('Excursion 1000Hz')
xlabel('Seconds')
ylabel('Milimeter')

% Values from plot
fi1 = 0.01143e-3;
fi2 = 0.0007973e-3;
tau1 = 1.035; %[s]      %4.474e4; Samples
tau2 = 1.105; %[s]
m = 2.6657e-04;
l = 0.0267;

q=15; %periods
tau_d_m= tau2-tau1; %[s]
r = (2 * m)/tau_d_m * log(fi1/fi2) %[Ns/m]

r_a = r/l %[Ns/m^2] 