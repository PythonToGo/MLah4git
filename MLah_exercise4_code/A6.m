clc; close all; clear;

% Parameters
fs = 48000;         
T = 3;              % Total signal duration in seconds
t_start = 0.5;      % Signal start time
t_stop = 2.8;       % Signal stop time
N_freqs = 179;      % Number of frequencies between f_start and f_stop
f_start = 55;       % Start frequency in Hz
f_stop = 22000;     % Stop frequency in Hz
block_size = 2^12;  % Block size (FFT window length)
ampl = 1;           % Amplitude of the multisine signal

%  Initialize audio interface and device reader
deviceReader = audioDeviceReader;
deviceReader.Driver = 'CoreAudio';

BufferSize = 1024;
saving = false;
aPR = audioPlayerRecorder('Device', 'Aggregate Device', ...
    'SampleRate', fs,...
    'BitDepth', '16-bit integer',...
    'SupportVariableSize', true,...
    'BufferSize', BufferSize, ...
    'PlayerChannelMapping', [1,2], ...
    'RecorderChannelMapping', [1,2]);

% Generate multisine signal
[t, sig, ms_indices] = generate_multisine(N_freqs, f_start, f_stop, fs, block_size, T);
sig = sig * ampl;

% Apply Hanning window
t_ramp = 200e-3;    % Ramp-up and ramp-down duration
n_ramp = floor(t_ramp * fs);
hann = hanning(2 * n_ramp, 'periodic');
sig(1:n_ramp) = sig(1:n_ramp) .* hann(1:n_ramp);
sig(end-n_ramp+1:end) = sig(end-n_ramp+1:end) .* hann(n_ramp+1:end);

%% Play and record

[playData, recData, N_underrun, N_overrun] = play_rec(aPR, sig);
release(aPR);

t_recData = (0:size(recData, 1)-1)/fs;


% Cut out the data for analysis
n_start = floor(t_start * fs);
n_stop = floor(t_stop * fs);
rec = recData(n_start:n_stop-1, :);
t_rec = (0:length(rec)-1) / fs;

% Calculate the FFT of the input and output signals
fft_rec = fft(rec) / block_size;
fft_rec = 2 * fft_rec(2:floor(block_size/2)+1, :);  % Single-sided spectrum without DC
fft_freqs = fs * (1:block_size/2)' / block_size;   % Frequencies of the spectrum

% Select the frequency bins where energy was provided
ms_fft_freqs = fft_freqs(ms_indices, :);
ms_fft_rec = fft_rec(ms_indices, :);

% Calculate the crosstalk CT12 at each frequency bin
CT12 = -20 * log10(abs(ms_fft_rec(:, 1)) ./ abs(ms_fft_rec(:, 2)));

fprintf('Crosstalk (f = 8 kHz): %.2f dB\n', CT12);

figure(1);
plot(CT12);
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Crosstalk CT12 (dB)','FontSize',14);
title('Crosstalk CT12 for R_{vict} = 1kohm','FontSize',16);
grid on;

%% Calculate transfer function
H=ms_fft_rec(:, 2) ./ ms_fft_rec(:, 1);
f = logspace(0, log10(24000), 1000);
w = 2*pi*f; % Converting frequency -> angular frequency
s = 1i*w;  

% Trnasfer function
H = @(s) 100000000 ./ (s.^2 + 30000.*s + 100000000);

% magnitude estimation
mag = abs(H(1j*w));
magDB = 20 * log10(mag);

% phase estimation
phaseDeg = unwrap(rad2deg(angle(H(1j*w))));

% plot
figure(2)
semilogx(w/(2*pi),magDB);
xlabel('frequency [Hz]','FontSize', 14);
ylabel('magniude & Crosstalk CT12 [dB]','FontSize', 14);
grid on;

hold on;

semilogx(ms_fft_freqs, CT12, 'r*');
title('Crosstalk dominating measurement','FontSize',16);
grid on;

