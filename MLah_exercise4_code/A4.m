%clc; clear; close all;

% generatea sine wave
Fs = 48000; % sample freuq
T = 1; % duration 
f = 8000; % sine freq
t = (0:1/Fs:T-1/Fs); % time vector
A = 0.7; % Amplitude
signal = A*sin(2*pi*f*t); % total wave

% output signal on channels
player = audioplayer([signal, signal], Fs); %????????????????????

% recoding signal in input CH1
recorder = audiorecorder(Fs, 16, 1); % to CH1 
recordblocking(recorder, T);
data_in = getaudiodata(recorder);

% crosstalk calculate 
A_out = max(signal); % Max ampl output open output
A_in = max(data_in(:,1)); % Max ampl input signal
CTo1 = -20*log10(A_out/A_in); % Crosstalk in dB in scripts

fprintf('Crosstalk (f = 8 kHz): %.2f dB\n', CTo1);
