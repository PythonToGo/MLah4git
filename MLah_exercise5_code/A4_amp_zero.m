clc, close all, clear all

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

% Parameters
T = 10;            % Total signal duration in (s)
t_start = 0.4;     % Adjust based on your system
t_stop = 9.6;      % Adjust based on your system
block_size = 2^12; % Block size (FFT window length)

amp = 0;       % 1 mV peak amplitude

N_freqs = 179;     % Number of frequencies
f_start = 55;      % Start frequency (Hz)
f_stop = 22000;    % Stop frequency (Hz)

f = 375;            % set signal frequency in Hz
N = floor(T*fs);    % calculate number of samples
t = (0:N-1)'/fs;    % generate time vector in s

%%
% zero-mean white noise       
white_noise = (rand(N, 1) * amp - amp / 2);


% plus a sig with amp and frequency to the noise
sinusoidal_signal = 0.05*sin(2*pi*f*t);
white_noise = white_noise + sinusoidal_signal; %% 그냥 플러스?

% hanning window
t_ramp = 0.01;    
n_ramp = floor(t_ramp * fs);
hann = hanning(2 * n_ramp, 'periodic');

% white noise
white_noise(1:n_ramp) = white_noise(1:n_ramp) .* hann(1:n_ramp);
white_noise(end-n_ramp+1:end) = white_noise(end-n_ramp+1:end) .* flip(hann(n_ramp+1:end));

% extend white noise to two channels and ch2 mute
white_noise(:,2) = 0;


%% Generate a multi-sine output signal
[t, white_noise, ms_indices] = generate_multisine(N_freqs, f_start, f_stop, fs, block_size, T);
sig = white_noise * amp;  % scale signal

% ramp-up, ramp-down with hanning window
t_ramp = 100e-3;
n_ramp = floor(t_ramp*fs);
hann = hanning(2*n_ramp, 'periodic');
sig(1:n_ramp) = sig(1:n_ramp) .* hann(1:n_ramp);
sig(end-n_ramp+1:end) = sig(end-n_ramp+1:end) .* hann(n_ramp+1:end);

%% Play and record
% Play and record (the recorded data will be the system noise u2N)
[playData, recData, N_underrun, N_overrun] = play_rec(aPR, sig*amp);
release(aPR);

% Define t_recData here, after the modification of recData
t_recData = (0:size(recData, 1)-1)'/fs;


%% Process measurement data 
n_start = floor(t_start*fs);
n_stop = floor(t_stop*fs);
num_avg = floor((n_stop-n_start+1) /block_size); % number of blocks in data
n_stop = n_start + num_avg*block_size;          % align end to block border

rec = recData(n_start:n_stop-1,:);  % cut out data for analysis

%% Plot the recorded noise (u2N)
figure;
plot(t_recData, recData(:, 2));
title('System Noise u2N(t)','FontSize',16);
xlabel('Time (s)','FontSize',14);
ylabel('Amplitude','FontSize',14);

% Noise floor with SNR
noise_floor = std(rec(:, 2));
% Output the noise floor
disp(['Noise floor in time domain for the total output (u2N(t)) is ', num2str(noise_floor), ' V'])

%% noise floor from dB
% Average in time domain
noise = recData(n_start:n_stop-1,2);  % cut out data for analysis
rec_avg=mean(reshape(rec,block_size,num_avg,2),2);
noise_avg=mean(reshape(noise,block_size,num_avg),2);

% FFT of input and output
% fft of the input
fft_rec = fft(rec_avg)/block_size;
fft_rec = 2*fft_rec(2:floor(block_size/2)+1,:); 

fft_freqs = fs*(1:block_size/2)'/block_size; % frequencies of spectrum

% fft of noise
fft_noise = fft(noise_avg)/block_size;
fft_noise = 2*fft_noise(2:floor(block_size/2)+1); % single sided spectrum without DC

% Calculate the magnitude of the noise floor in dB
noise_floor = 20*log10(abs(fft_noise));


%% Plot noise floor

figure(442);

semilogx(fft_freqs, abs(fft_rec(:,2)), 'LineWidth',3);
title('Average Spectrum','FontSize',16);


hold on; 
semilogx(fft_freqs, noise_floor);
title('Noise Floor','FontSize',16);
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Magnitude (dB)','FontSize',14);
legend('Average Spectrum','Noise Floor','Fontsize', 14);
hold off;


