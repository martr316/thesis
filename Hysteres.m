% Hysteresis

 load('DCstep_5V_Jadetär5V');
% figure
% plot(MeasureddataExcersion.Data)
% title('5volt')
% 
% load('DCstep_10V_Jadetär10V');
% figure
% plot(MeasureddataExcersion.Data)
% title('10volt')
% 
% load('DCstep_15V_Jadetär15V');
% figure
% plot(MeasureddataExcersion.Data)
% title('15volt')
% 
% load('DCstep_20V_Jadetär20V');
% figure
% plot(MeasureddataExcersion.Data)
% title('20volt')
% 
% load('DCstep_25V_Jadetär25V');
% figure
% plot(MeasureddataExcersion.Data)
% title('25volt')
% 
% load('DCstep_30V_Jadetär30V');
% figure
% plot(MeasureddataExcersion.Data)
% title('30volt')

%% Hysteresis 5 volt
load('DCstep_5V_Jadetär5V');
up_sample5=1:(6.575e4-6.341e4);
down_sample5=1:(2.222e5-2.209e5);
data_up5=0.008224+MeasureddataExcersion.Data(6.341e4:6.575e4); %calibrate to start from zero
data_down5=0.006553+MeasureddataExcersion.Data(2.209e5:2.222e5); %calibrate to start from zero

time_total5 = [down_sample5 fliplr(down_sample5)]; %shortest time vector;
sampled_data_up5 = data_up5(1:length(data_up5)/length(data_down5):end);
data_total5 = [sampled_data_up5(1:(end-2)); data_down5]; %shorten with two to fit vectors

figure(50)
plot(time_total5,data_total5)
title('5 volt')
ylabel('Excusion [mm]')

%% Hysteresis 10 volt
load('DCstep_10V_Jadetär10V');
up_sample10=1:(5.727e4-5.546e4);
down_sample10=1:(1.853e5-1.828e5);
data_up10=0.0155+MeasureddataExcersion.Data(5.546e4:5.727e4); %calibrate to start from zero
data_down10=0.01305+MeasureddataExcersion.Data(1.828e5:1.853e5); %calibrate to start from zero

time_total10 = [down_sample10 fliplr(down_sample10)]; %shortest time vector;
sampled_data_up10 = data_up10(1:length(data_up10)/length(data_down10):end);
data_total10 = [sampled_data_up10(1:(end-1)); data_down10]; %shorten with two to fit vectors

figure(51)
plot(time_total10,data_total10)
title('10 volt')
ylabel('Excusion [mm]')

%% Hysteresis 15 volt
load('DCstep_15V_Jadetär15V');
up_sample15=1:(1.07e5-1.053e5);
down_sample15=1:(2.546e5-2.51e5);
data_up15=0.02793+MeasureddataExcersion.Data(1.053e5:1.07e5); %calibrate to start from zero
data_down15=0.02632+MeasureddataExcersion.Data(2.51e5:2.546e5); %calibrate to start from zero

time_total15 = [down_sample15 fliplr(down_sample15)]; %shortest time vector;
sampled_data_up15 = data_up15(1:length(data_up15)/length(data_down15):end);
data_total15 = [sampled_data_up15; data_down15];

figure(52)
plot(time_total15,data_total15)
title('15 volt')
ylabel('Excusion [mm]')

%% plot for all hysteresis
figure(59)
plot(time_total5,data_total5, time_total10,data_total10, time_total15,data_total15)

%%


diff_ex5 = 0.009867+0.01611;
diff_ex10 = 0.01588 + 0.03735;
diff_ex15 = 0.02833 +0.05491;
diff_ex20 = 0.06596 + 0.04973;
diff_ex25 = 0.05543 + 0.09415;
diff_ex30 = 0.06847 + 0.117;


voltage = [0 5 10 15 20 25 30];
excursion = [0 diff_ex5 diff_ex10 diff_ex15 diff_ex20 diff_ex25 diff_ex30];

plot(voltage,excursion);
xlabel('Voltage [v]')
ylabel('Excursion [mm]')