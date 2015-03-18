samplerate=44100;
secondes = 5;
t = 0:1:samplerate*3;  % Time Samples
f = 4300;     % Input Signal Frequency
fs = 44100; % Sampling Frequency


stepvector1 = [sin(2*pi*f/fs*t) zeros(1,samplerate)];
%stepvector2 = [2*ones(1,samplerate*5) zeros(1,samplerate*10)];
%stepvector5 = [5*ones(1,samplerate*5) zeros(1,samplerate*10)];
%stepvector10 = [10*ones(1,samplerate*5) zeros(1,samplerate*10)];
plot(stepvector1);

%wavwrite(stepvector1,samplerate,'test2')

audiowrite( 'Step_1_0.wav',stepvector1, samplerate);
%audiowrite( 'Step_2_0.wav',stepvector2, samplerate);
%audiowrite( 'Step_5_0.wav',stepvector5, samplerate);
%audiowrite( 'Step_10_0.wav',stepvector10, samplerate);
