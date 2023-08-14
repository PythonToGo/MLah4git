%% __INIT__
clc             % clear command line
clear all       % clear all variables, functions, ... to avoid random conditions
close all       % close all open windows

% Setup Audio Interface
deviceReader = audioDeviceReader;
deviceReader.Driver = 'CoreAudio';
%devices = getAudioDevices(deviceReader);

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
T = 5;             % Total signal duration in s
t_start = 0.4;      % signal analysis start: remove fading and transient part
t_stop=4.8;           % signal analysis end: remove the fading part

f = 100;            % set signal frequency in Hz
N = floor(T*fs);    % calculate number of samples
t = (0:N-1)'/fs;    % generate time vector in s
amp = 0.5;          % set peak amplitude of output re full scale <=1

% Generate impulse signal
impulse = zeros(N, 1);
impulse(round(N/2)) = 1; % Set the amplitude to 1 at the middle point

% sine wave signal
sig = amp*sin(2*pi*f*t);            % calculate sine signal
sig(:,2) = amp*sin(2*pi*f*t);       % set signal for second channel

% sine wave signal
white_noise = rand(N, 1) * amp - amp / 2;          
%white_noise(:,1) = rand(N, 1) * amp - amp / 2;% set signal for second channel


% Apply Hanning window
t_ramp = 0.1;    % Ramp-up and ramp-down duration
n_ramp = floor(t_ramp * fs);
hann = hanning(2 * n_ramp, 'periodic');
white_noise(1:n_ramp) = white_noise(1:n_ramp) .* hann(1:n_ramp);
white_noise(end-n_ramp+1:end) = white_noise(end-n_ramp+1:end) .* flip(hann(n_ramp+1:end));

% extend white oise to two channels and ch2 mute
white_noise(:,2) = 0;


%% Play and Record signal
fprintf('Start recording  -- read output voltage with voltmeter (RMS)\n...\n')
[playData, recData] = play_rec(aPR, white_noise);
fprintf('End recording\n')
t_recData = (0:length(recData)-1)/fs;

release(aPR);
%rms_amp=input('RMS amplitude in V_rms:');    % input RMS value
rms_amp = 0.5;
fprintf('this coresponds to a peak amplitude of %.4f V\n',rms_amp*sqrt(2));

% Cut out the data for analysis
n_start = floor(t_start * fs);
n_stop = floor(t_stop * fs);
rec = recData(n_start:n_stop-1, :);
t_rec = (0:length(rec)-1) / fs;


%% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%Assignment1 (a)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
% Plot u1(t)
subplot(211); 
plot(t_recData, recData(:,1)); 
title('Signal u1(t)','FontSize',18);
xlabel('Time (s)','FontSize',16);
ylabel('Amplitude','FontSize',16);
% Plot u2(t)
subplot(212); 
plot(t_recData, recData(:,2)); 
title('Signal u2(t)','FontSize',18);
xlabel('Time (s)','FontSize',16);
ylabel('Amplitude','FontSize',16);

%%
% Assuming rec(:,1) contains u1(t) and rec(:,2) contains u2(t)
mean_u1 = mean(rec(:,1));
std_u1 = std(rec(:,1));

mean_u2 = mean(rec(:,2));
std_u2 = std(rec(:,2));

fprintf('For u1(t), mean: %f, standard deviation: %f\n', mean_u1, std_u1);
fprintf('For u2(t), mean: %f, standard deviation: %f\n', mean_u2, std_u2);
