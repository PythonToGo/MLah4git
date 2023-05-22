%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab Funcion for Measurement Laboratory at Home
% The script generates a sine signal, ouputs it and reads the inputs back
% Needs the subroutine play_rec.m
% aPR is the MATLAB audioPlayerRecorder
% You might have to adjust at least the paramters t_start and t_stop such 
% that the code works for your system.
% Change the amplitude amp at your needs.
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
% NOTE: the audio interface needs an ASIO4all driver!!!
% by Lukas Driendl, Bernhard Gleich, Ali Saeedi & Werner Hemmert, TUM, 
% 20. Sept. 2020 - 14. Oct. 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc             % clear command window
clear all       % clear all variables, functions, ... to avoid random conditions
close all       % close all open windows

% Setup Audio Interface

deviceReader.Driver = 'CoreAudio';
deviceReader = audioDeviceReader;

%%devices = getAudioDevices(deviceReader);

% Define BufferSize and Sampling Frequency (has to match ASIO Settings)
BufferSize = 1024;      % ASIO Buffer Size in Samples
fs = 48000;             % Sampling frequency in Hz

% Create Audio Player Recorder object
aPR = audioPlayerRecorder('Device', 'Aggregate Device',...
    'SampleRate', fs,...
    'BitDepth', '16-bit integer',...
    'SupportVariableSize', true,...
    'BufferSize', BufferSize, ...
    'PlayerChannelMapping', [1,2], ...
    'RecorderChannelMapping', [1,2]);

%% Parameters
T = 2;              % Total signal duration in s
t_start = 0.4;      % signal analysis start: remove fading and transient part
t_stop=1.4;         % signal analysis end: remove fading part

f = 100;            % set signal frequency in Hz
N = floor(T*fs);    % calculate number of samples
t = (0:N-1)'/fs;    % generate time vector in s
amp = 1;            % set peak amplitude of output re full scale <=1

sig = amp*sin(2*pi*f*t);            % calculate sine signal
sig(:,2) = 1*amp*sin(2*pi*f*t);     % set signal for second channel to same signal

% ramp-up, ramp-down signal with hanning window
t_ramp = 100e-3;                    % ramp-up and down time in s
n_ramp = floor(t_ramp*fs);          % calculate number of samples for ramp
hann = hanning(2*n_ramp, 'periodic');     % ramp is a hanning window shape
hann(:,2) = hanning(2*n_ramp, 'periodic');
sig(1:n_ramp,:) = sig(1:n_ramp,:) .* hann(1:n_ramp,:); % apply ramp on
sig(end-n_ramp+1:end,:) = sig(end-n_ramp+1:end,:) .* hann(n_ramp+1:end,:);

% Play and Record signal
fprintf('Start recording...\n')
[playData, recData] = play_rec(aPR, sig);
fprintf('End recording\n')
t_recData = (0:length(recData)-1)/fs;

% Release Audio Player Recorder Object
release(aPR);

%% Plot input signal
figure
subplot(2,2,1)
plot(t_recData, recData(:, 1), 'k');
ylabel('CH1 input /a.u.');
title('Raw signal');
subplot(2,2,3)
plot(t_recData, recData(:, 2), 'r');
xlabel('Time/s');
ylabel('CH2 input /a.u.');

% zoom into signal from t_start on 5 periods
n_start = floor(t_start*fs);
n_stop = n_start+floor(5*fs/f);   % 5 periods

subplot(2,2,2)
plot(t_recData(n_start:n_stop), recData(n_start:n_stop, 1), 'k');
ylabel('CH1 input /a.u.');
title('Raw signal (zoom in)');
subplot(2,2,4)
plot(t_recData(n_start:n_stop), recData(n_start:n_stop, 2), 'r');
xlabel('Time/s');
ylabel('CH2 input /a.u.');

% print('check_uai', '-dpng', '-r300');     % uncomment to save figure as .png
