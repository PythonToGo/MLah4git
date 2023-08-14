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
T = 10;             % Total signal duration in s
t_start = 0.4;      % signal analysis start: remove fading and transient part
t_stop = 9.6;       % signal analysis end: remove the fading part

f = 375;            % set signal frequency in Hz
N = floor(T*fs);    % calculate number of samples
t = (0:N-1)'/fs;    % generate time vector in s

% zero-mean white noise
amp = 0.5;          
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

% sinus sigal
%{
sinusoidal_signal(1:n_ramp) = sinusoidal_signal(1:n_ramp) .* hann(1:n_ramp);
sinusoidal_signal(end-n_ramp+1:end) = sinusoidal_signal(end-n_ramp+1:end) .* flip(hann(n_ramp+1:end));
%}

% extend white noise to two channels and ch2 mute
white_noise(:,2) = 0;
%sinusoidal_signal(:,2) = 0;

%% Play and Record signal
fprintf('Start recording  -- read output voltage with voltmeter (RMS)\n...\n')
[playData, recData] = play_rec(aPR, white_noise);
fprintf('End recording\n')
t_recData = (0:length(recData)-1)/fs;

release(aPR);
%rms_amp=input('RMS amplitude in V_rms:');    % input RMS value
rms_amp = 0.5;
fprintf('this coresponds to a peak amplitude of %.4f V\n',rms_amp*sqrt(2));

% Extract u1(t) and u2(t) from recorded data
u1 = recData(:, 1);
u2 = recData(:, 2);

% plot recorded Data
%{
figure(1);
plot(t_recData, recData);
%}

%% Plot u1 u2
%{
% Plot u1(t)
figure(2);
subplot(2,1,1);
plot(t_recData, u1);
title('u1(t)','Fontsize', 16);
xlabel('Time [s]', 'Fontsize', 14);
ylabel('Amplitude', 'Fontsize', 14);

% Plot u2(t)
subplot(2,1,2);
plot(t_recData, u2);
title('u2(t)','Fontsize', 16);
xlabel('Time [s]','FontSize',14);
ylabel('Amplitude','FontSize',14);
%}

%% Divide recSignal to make average

% Define block size
block_size = 2^12;

% Extract part of the signal for analysis
n_start = floor(t_start*fs);
n_stop = floor(t_stop*fs);
num_avg = floor((n_stop-n_start+1) /block_size); % number of blocks in data
n_stop = n_start + num_avg*block_size;          % align end to block border

rec = recData(n_start:n_stop-1,:);  % cut out data for analysis

% Average in time domain
rec_avg = mean(reshape(rec,block_size,num_avg,2),2);
t_avg = (0:block_size-1)'/fs;  % Time vector for the averaged signal

%% Plot averaged signals
%{
% Plot the averaged signal for u1
figure(3);
subplot(2, 1, 1);
plot(t_avg, rec_avg(:, 1));
title('Averaged Signal u1_{avg}(t)','Fontsize', 16);
xlabel('Time [s]', 'Fontsize', 14);
ylabel('Amplitude', 'Fontsize', 14);

% Plot the averaged signal for u2
subplot(2, 1, 2);
plot(t_avg, rec_avg(:, 2));
title('Averaged Signal u2_{avg}(t)','Fontsize', 16);
xlabel('Time [s]', 'Fontsize', 14);
ylabel('Amplitude', 'Fontsize', 14);
%}

%% Calculate Mean and Standard Deviation

% For u1_avg(t) : mean std
mean_u1_avg = mean(rec_avg(:, 1));
std_u1_avg = std(rec_avg(:, 1));

%{
fprintf('The mean of u1_avg is: %.4f\n', mean_u1_avg);
fprintf('The standard deviation of u1_avg is: %.4f\n', std_u1_avg);
%}

% For u2_avg(t) : mean std
mean_u2_avg = mean(rec_avg(:, 2));
std_u2_avg = std(rec_avg(:, 2));

%{
fprintf('The mean of u2_avg is: %.4f\n', mean_u2_avg);
fprintf('The standard deviation of u2_avg is: %.4f\n', std_u2_avg);
%}

%% FFT U2(f)
% FFT for each block
u2 = u2(n_start:n_stop-1);
u2_blocks = reshape(u2, [block_size, num_avg]);
U2_f_blocks = fft(u2_blocks);

%{
Now, U2_f_blocks == a 2d matrics, (each column is the FFT)
For magnitude: abs(U2_f_blocks)
For phase: angle(U2_f_blocks)
%}

% frequency vector
f_vec = (0:block_size-1)*fs/block_size;

% Plot %%%%%%%%%%%%%%%%% A3(a) %%%%%%%%%%%%%%%%%%%%
figure(31);
plot(f_vec(1:block_size/2), abs(U2_f_blocks(1:block_size/2, 1)));
title('Magnitude Spectrum of U2(f)','FontSize',16);
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Magnitude','FontSize',14);

%% %%%%%%%%%%%%%%% A3(b) %%%%%%%%%%%%%%%%%%%%%%
% FFT and FFT_avg
U2_blocks_fft = fft(u2_blocks);
U2_blocks_fft_avg = mean(U2_blocks_fft, 2);

% magnitude
U2_blocks_fft_avg_mag = abs(U2_blocks_fft_avg);

% frequency vector
f_vec = (0:block_size-1)*fs/block_size;

% Plot
figure(32);
plot(f_vec, 20*log10(U2_blocks_fft_avg_mag));
title('Magnitude of Averaged Complex Spectra','FontSize',16);
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Magnitude (dB)','FontSize',14);

%% %%%%%%%%%%%%%% A3(c) %%%%%%%%%%%%%%%%%%%%%%%%%
% magnitude of FFT
U2_blocks_fft_mag = abs(U2_blocks_fft);

% Average!: magnitude
U2_blocks_fft_mag_avg = mean(U2_blocks_fft_mag, 2);

% Plot
figure(33);
plot(f_vec, 20*log10(U2_blocks_fft_mag_avg));
title('Average of Magnitude of Spectra with ','FontSize',16);
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Magnitude (dB)','FontSize',14);

%% %%%%%%%%%%%%%%%%% A3(d) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
block_sizes = [2^9, 2^10, 2^11];
colors = ['r', 'g', 'b'];
figure(34);

for i = 1:length(block_sizes)
    block_size = block_sizes(i);
    num_avg = floor(length(u2) / block_size);
    
    % reshape, compute FFT
    u2_blocks = reshape(u2(1:block_size*num_avg), [block_size, num_avg]);
    U2_blocks_fft = fft(u2_blocks);
    U2_blocks_fft_mag = abs(U2_blocks_fft);
    
    % Average!  magnitudes
    U2_blocks_fft_mag_avg = mean(U2_blocks_fft_mag, 2);
    
    % frequency vector
    f_vec = linspace(0, fs, block_size);
    
    % Plot 
    semilogx(f_vec, 20*log10(U2_blocks_fft_mag_avg), colors(i));
    hold on;
end

title('Average of Magnitude of Spectra for Every Block Sizes','FontSize',16);
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Magnitude (dB)','FontSize',14);
legend('Block size = 2^9', 'Block size = 2^{10}', 'Block size = 2^{11}','Fontsize', 14);
hold off;
