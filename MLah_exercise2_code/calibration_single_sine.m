%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab Funcion for Measurement Laboratory at Home
% The script generates a sine signal, ouputs it and reads the inputs back
% Needs the subroutine play_rec.m
% aPR is the MATLAB audioPlayerRecorder
% You will have to adjust at least the paramters t_start and t_stop such 
% that the code works for your system.
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
% NOTE: the audio interface needs an ASIO4all driver!!!
% by Lukas Driendl, Ali Saeedi & Werner Hemmert, TUM, 20. Sept. 2020 - 14. Oct. 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc             % clear command line
clear all       % clear all variables, functions, ... to avoid random conditions
close all       % close all open windows

% Setup Audio Interface

deviceReader.Driver = 'CoreAudio';
deviceReader = audioDeviceReader;
devices = getAudioDevices(deviceReader);

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
T = 5;              % Total signal duration in s
t_start = 0.4;      % signal analysis start: remove fading and transient part
t_stop = 4.6;         % signal analysis end: remove the fading part

f = 100;            % set signal frequency in Hz
N = floor(T*fs);    % calculate number of samples
t = (0:N-1)'/fs;    % generate time vector in s
amp = 1.0;          % set peak amplitude of output re full scale <=1

white_noise = rand(N, 1) * 0.5;     % calculate rand white noise
white_noise(:,2) = rand(N, 1) * 0.5;     % set signal for second channel

% zero mean white noise signal
T = 0.5;
amplitude = 0.5;

round_trip_log = [];

for i = 1:10
    % Generate white noise signal
    
    numSamples = T * fs; % Number of samples
    new_white_noise_signal = amplitude * randn(numSamples, 1); % Generate white noise signal with zero mean
    new_white_noise_signal(:,2) = amplitude * randn(numSamples, 1);
    
    sig = amp*sin(2*pi*f*t);            % calculate sine signal
    sig(:,2) = amp*sin(2*pi*f*t);       % set signal for second channel
    
    %% Generate impulse signal vector
    
    % Set the parameters
    duration = 0.25;  % Duration of the signal in seconds
    
    % Create a zero vector
    numSamples = round(duration * fs);
    impulse = zeros(numSamples, 1);
    
    % Set the sample at 20 ms to 1
    sampleIndex = round(0.02 * fs); % 20 ms
    impulse(sampleIndex) = 1;
    
    
    
    %%
    % ramp-up, ramp-down signal with hanning window
    t_ramp = 100e-3;                    % ramp-up and down time in s
    n_ramp = floor(t_ramp*fs);          % calculate number of samples for ramp
    hann = hanning(2*n_ramp, 'periodic');     % ramp is a hanning window shape
    hann(:,2) = hanning(2*n_ramp, 'periodic');
    sig(1:n_ramp,:) = sig(1:n_ramp,:) .* hann(1:n_ramp,:); % apply ramp on
    sig(end-n_ramp+1:end,:) = sig(end-n_ramp+1:end,:) .* hann(n_ramp+1:end,:);
    
    
    new_white_noise_signal(1:n_ramp,:) = new_white_noise_signal(1:n_ramp,:) .* hann(1:n_ramp,:); % apply ramp on
    new_white_noise_signal(end-n_ramp+1:end,:) = new_white_noise_signal(end-n_ramp+1:end,:) .* hann(n_ramp+1:end,:);
    
    
    % set input signal here
    input_signal = new_white_noise_signal;
    
    
    % Play and Record signal



    
    fprintf('Start recording  -- read output voltage with voltmeter (RMS)\n...\n')
    [playData, recData] = play_rec(aPR, input_signal);
    fprintf('End recording\n')
    t_recData = (0:length(recData)-1)/fs;
    
    release(aPR);
    %rms_amp=input('RMS amplitude in V_rms:');   % input RMS value
    rms_amp=0.7; 
    fprintf('this coresponds to a peak amplitude of %.4f V\n',rms_amp*sqrt(2));

%% Plot white noise
%{
figure
subplot(2,1,1)
plot(t_recData, playData,'k');  % scale peak amplitude
ylim([-amp-amp/10, amp+amp/10])
ylabel('CH1 Amplitude in V');
xlabel('Time in s');
title('Output signal of Playdata (scaled in Volt)');

subplot(2,1,2)
plot(t_recData, recData,'r');  % assume same gain as CH1
ylim([-amp-amp/10, amp+amp/10])
ylabel('CH2 Amplitude in V');
xlabel('Time in s');
title('Output signal of Recdata (scaled in Volt)');


%% Plot scaled output signal
figure
subplot(2,1,1)
plot(t, sig(:,1)*rms_amp*sqrt(2)/amp,'k');  % scale peak amplitude
ylim([-amp-amp/10, amp+amp/10])
ylabel('CH1 Amplitude in V');
title('Output signal (scaled in Volt)');
subplot(2,1,2)
plot(t, sig(:,2)*rms_amp*sqrt(2)/amp,'r');  % assume same gain as CH1
ylim([-amp-amp/10, amp+amp/10])
ylabel('CH2 Amplitude in V');
xlabel('Time in s');


%% Plot input signal (unscaled, raw signal)

figure
subplot(2,2,1)
plot(t_recData, recData(:, 1), 'k');
ylabel('CH1 Amplitude re full scale');
title('Measured raw signal');
subplot(2,2,3)
plot(t_recData, recData(:, 2), 'r');
xlabel('Time/s');
ylabel('CH2 Amplitude re full scale');

% zoom into signal from t_start on 5 periods
n_start = floor(t_start*fs);
n_stop = n_start+floor(5*fs/f);   % 5 periods

subplot(2,2,2)
plot(t_recData(n_start:n_stop), recData(n_start:n_stop, 1), 'k');
ylabel('CH1 Amplitude re full scale');
title('Measured raw signal (zoom in)');
subplot(2,2,4)
plot(t_recData(n_start:n_stop), recData(n_start:n_stop, 2), 'r');
xlabel('Time/s');
ylabel('CH2 Amplitude re full scale');


%% Generate impulse signal vector

% Set the parameters
duration = 0.25;  % Duration of the signal in seconds

% Create a zero vector
numSamples = round(duration * fs);
impulse = zeros(numSamples, 1);

% Set the sample at 20 ms to 1
sampleIndex = round(0.02 * fs); % 20 ms
impulse(sampleIndex) = 1;
% set input signal here
input_signal = impulse;

% Plotting the signals from recData and playData
figure;
subplot(2, 1, 1);
plot(playData);
title('Played Signal');
xlabel('Sample Index');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(recData);
title('Recorded Signal');
xlabel('Sample Index');
ylabel('Amplitude');

% Zooming into the played and recorded impulses
zoomStart = round(0.02 * fs); % Starting index for zooming (20 ms)
zoomEnd = round(0.03 * fs); % Ending index for zooming (30 ms)

figure;
subplot(2, 1, 1);
plot(playData(zoomStart:zoomEnd));
title('Zoomed Played Signal');
xlabel('Sample Index');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(recData(zoomStart:zoomEnd));
title('Zoomed Recorded Signal');
xlabel('Sample Index');
ylabel('Amplitude');

%{
% Find the index of the maximum value in the recorded signal
[~, maxIndex] = max(recData);

% Calculate the round-trip delay in number of samples
roundTripDelaySamples = maxIndex - sampleIndex;

% Calculate the round-trip delay in milliseconds
roundTripDelayMilliseconds = roundTripDelaySamples / fs * 1000;

fprintf('the round-trip delay is:  %.4f  \n',roundTripDelayMilliseconds)
%}
%}

%% Plot the cross-correlation
% zero mean white noise signal
    T = 0.5;
    amplitude = 0.5;
    
    % Generate white noise signal
    
    numSamples = T * fs; % Number of samples
    signal = amplitude * randn(numSamples, 1); % Generate white noise signal with zero mean
    
    % Compute the cross-correlation between the output and input signals
    [correlation, lags] = xcorr(recData(:,1), input_signal(:,1));
    
    % Plot the cross-correlation
    plot(lags, correlation);
    xlabel('Lag (samples)');
    ylabel('Cross-Correlation');
    title('Cross-Correlation between Output and Input Signals');
    
    % Find the index of the maximum value in the recorded signal
    [~, maxIndex] = max(recData(:,1));
    
    % Calculate the round-trip delay in number of samples
    roundTripDelaySamples = maxIndex - sampleIndex;
    
    % Calculate the round-trip delay in milliseconds
    roundTripDelayMilliseconds = roundTripDelaySamples / fs * 1000;
    
    fprintf('the round-trip delay is:  %.4f  \n',roundTripDelayMilliseconds)

    round_trip_log = [round_trip_log, roundTripDelayMilliseconds];

end

% Display the round-trip delays
disp("Round-Trip Delays:")
disp(round_trip_log)

% Calculate the jitter
jitter = max(round_trip_log) - min(round_trip_log);

% Display the jitter
disp("Jitter:")
disp(jitter)

%% Analyze, plot selected input signal part

% remove the transient part
n_start = floor(t_start*fs);
n_stop = floor(t_stop*fs);

rec = recData(n_start:n_stop-1, :);
t_rec = (0:length(rec)-1)/fs;

figure;
subplot(2,1,1)
plot(t_rec, rec(:, 1), 'k');
ylabel('CH1 Amplitude re full scale');
title('Measured signal for analysis(unscaled)');
subplot(2,1,2)
plot(t_rec, rec(:, 2), 'r');
xlabel('Time/s');
ylabel('CH2 Amplitude re full scale');

%% FFT of input
N = length(rec);
ss_fft_freqs = fs*(0:N/2)/N;

[min_diff, f_index] = min(abs(ss_fft_freqs - f));

% fft of the input
fft_rec = fft(rec)/N;
ss_fft_rec = [fft_rec(1,:); 2*fft_rec(2:floor(N/2)+1,:)];

% extract amplitude
a_f = abs(ss_fft_rec(f_index,:));  % amplitude at f of CH1 and CH2

%% Calculate Output and Input sensitivity
% Output Sensitivity
sensitivity_out = amp ./ (rms_amp * sqrt(2));

% Input Sensitivity
sensitivity_in = rms_amp * sqrt(2) ./ a_f;

% Save calibration for later use
save('sensitivity.mat', 'sensitivity_in', 'sensitivity_out');  
fprintf('\n-----------------------------------------------------------------\n')
fprintf('The output sensitivity of your output channels is Ch1&2: %.4f (in full scale/V_p) \n', sensitivity_out)
fprintf('The input sensitivity of your input channels is CH1: %.4f  CH2: %.4f (in V_p/full scale)\n',sensitivity_in)
fprintf('-----------------------------------------------------------------\n')

%% Plot the Input Spectrum with correct scaling in Volts
rec_scaled = rec .* sensitivity_in;     % now always scale your inputs!!

% FFT of scaled input signal
fft_rec_scaled = fft(rec_scaled)/N;
ss_fft_rec_scaled = [fft_rec_scaled(1,:); 2*fft_rec_scaled(2:floor(N/2)+1,:)];

% Plot the spectra
figure
subplot(2,2,1)
stem(ss_fft_freqs, abs(ss_fft_rec_scaled(:, 1)), 'k');
ylabel('CH1 Magnitude in V');
title('Scaled amplitude spectrum (lin)');
set(gca,'xscal','log')
xticks([1e1, 1e2, 1e3, 1e4])
subplot(2,2,3)
stem(ss_fft_freqs, abs(ss_fft_rec_scaled(:, 2)), 'r');
xlabel('f/Hz');
ylabel('CH2 Magnitude in V');
xticks([1e1, 1e2, 1e3, 1e4])
set(gca,'xscal','log')

subplot(2,2,2)
loglog(ss_fft_freqs, abs(ss_fft_rec_scaled(:, 1)), 'k');
ylabel('CH1 Magnitude in V');
axis([10 fs/2 1e-7 3])
set(gca,'xscal','log')
xticks([1e1, 1e2, 1e3, 1e4])
title('Scaled amplitude spectrum (log)');
subplot(2,2,4)
loglog(ss_fft_freqs, abs(ss_fft_rec_scaled(:, 2)), 'r');
xlabel('f/Hz');
ylabel('CH2 Magnitude in V');
axis([10 fs/2 1e-7 4])
xticks([1e1, 1e2, 1e3, 1e4])
set(gca,'xscal','log')
