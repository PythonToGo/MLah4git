%clc; clear; close all;

% generatea sine wave
Fs = 48000; % sample freuq
T = 1; % duration 
f = 8000; % sine freq
t = (0:1/Fs:T-1/Fs); % time vector
A = 0.7; % Amplitude
signal = A*sin(2*pi*f*t); % total wave

% Output signal on channel 1 and input signal on channel 2
player = audioplayer([signal; zeros(1, length(t))], Fs);  
recorder = audiorecorder(Fs, 16, 2);                      % Output on channel 2
recordblocking(recorder, T);                             
data_in = getaudiodata(recorder);                        

% Calculate crosstalk CT12
A_out = max(signal);        % Max ampl output signal on channel 1
%A_out = max(signal);
A_in = max(data_in(:, 2));         % Max ampl input signal on channel 2
CT12 = -20*log10(A_out/A_in);      % Crosstsalkl

fprintf('Crosstalk (f = 8 kHz): %.2f dB\n', CT12);
