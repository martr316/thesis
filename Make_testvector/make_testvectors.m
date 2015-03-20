samplerate=44100;
secondes = 5;
t = 0:1:samplerate*3;  % Time Samples
f = 4300;     % Input Signal Frequency
fs = 44100; % Sampling Frequency


stepvector1 = [sin(2*pi*f/fs*t) zeros(1,samplerate)];
plot(stepvector1);
audiowrite( 'Step_1_0.wav',stepvector1, samplerate);

%% Create a musltisin205and1000Hz funtion
samplerate=44100;
seconds = 5;
t = 0:1:samplerate*seconds;  % Time Samples
f1 = 205;     % first resonance Frequency
f2 = 1000;    % low excursion Frequency
fs = 44100;   % Sampling Frequency

musltisin205and1000Hz = [sin(2*pi*f1/fs*t).*sin(2*pi*f2/fs*t) ];
plot(musltisin205and1000Hz);
audiowrite( 'musltisin205_1000Hz.wav',musltisin205and1000Hz, samplerate);

%% Create a musltisin1430and5000Hz funtion
samplerate=44100;
seconds = 5;
t = 0:1:samplerate*seconds;  % Time Samples
f1 = 1430;     % second resonance Frequency
f2 = 5000;   % low excursion Frequency

fs = 44100;   % Sampling Frequency

musltisin1430and5000Hz = [sin(2*pi*f1/fs*t).*sin(2*pi*f2/fs*t) ];
plot(musltisin1430and5000Hz);
audiowrite( 'musltisin1430_5000Hz.wav',musltisin1430and5000Hz, samplerate);

%% Create a 205 and 760 and 1430HZ funtion
samplerate=44100;
seconds = 5;
t = 0:1:samplerate*seconds;  % Time Samples
f1 = 205;     % first reconance Frequency
f2 = 760;     % first node freq
f3 = 1430;   % second resonance Frequency

fs = 44100;   % Sampling Frequency

musltisin205_760_1430Hz = [sin(2*pi*f1/fs*t).*sin(2*pi*f2/fs*t).*sin(2*pi*f3/fs*t)];
plot(musltisin205_760_1430Hz);
audiowrite( 'musltisin205_760_1430Hz.wav',musltisin205_760_1430Hz, samplerate);

%% Create a 205 and 760 and 1430HZ funtion
samplerate=44100;
seconds = 5;
t = 0:1:samplerate*seconds;  % Time Samples
f1 = 205;      % first reconance Frequency
f2 = 500;     
f3 = 760;      % first node freq
f4 = 1100;      
f5 = 1430;     % second resonance Frequency
f6 = 2500;     
f7 = 4315;     % third resonance Frequency
f8 = 6000;           

fs = 44100;   % Sampling Frequency

musltisin205_500_760_1100_1430_2500_4315_6000Hz = [sin(2*pi*f1/fs*t).*sin(2*pi*f2/fs*t).*sin(2*pi*f3/fs*t).*sin(2*pi*f4/fs*t).*sin(2*pi*f5/fs*t).*sin(2*pi*f6/fs*t).*sin(2*pi*f7/fs*t).*sin(2*pi*f8/fs*t)];
plot(musltisin205_500_760_1100_1430_2500_4315_6000Hz);
audiowrite( 'musltisin205_500_760_1100_1430_2500_4315_6000Hz.wav',musltisin205_500_760_1100_1430_2500_4315_6000Hz, samplerate);

