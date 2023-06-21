%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab Funcion for Measurement Laboratory at Home
% The script uses a multi-sine signal and writes and reads it
% needs the subroutines generate_multisine.m and play_rec.m
% aPR is the MATLAB audioPlayerRecorder
% You will have to adjust at least the paramters t_start and t_stop such 
% that the code works for your system.
% PS: If you use the multi sine method later in your career, do not 
% forget to cite our first publications:
% W. Hemmert, H.-P. Zenner, A.W. Gummer (1995): “Dreidimensionale Schwingungsmessung im Cortischen Organ der Säugetiercochlea”, FORTSCHRITTE DER AKUSTIK 21, p. 207-210.
% AW Gummer, W Hemmert, HP Zenner (1996): “Resonant tectorial membrane motion in the inner ear: its crucial role in frequency tuning”, Proceedings of the National Academy of Sciences 93 (16), p. 8727-8732.
% W. Hemmert (1997): “3D vibration measurements of the inner ear – Elucidation of frequency selectivity of the hearing organ”, PTB-MITTEILUNGEN 107 (1), p. 12-18. 
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
% NOTE: the audio interface needs an ASIO4all driver!!!
% by Ali Saeedi & Werner Hemmert, TUM, 20. Sept. 2020 - 6. Dec. 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, close all, clear all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  You have to set these values based on your individual system properties!
T = 3;            % Total signal duration in (s), more averages with longer signal
t_start = 0.5;      % signal start: remove also transient part
t_stop=2.8;         % signal end: remove the fading part
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Define parameters of your measurement
name='Test Transfer Function Measurement V0';  % name for output file
fs = 48000;             % sampling frequency, make sure value is correct!
N_freqs = 179;          % number of frequencies between f_start and f_stop 
f_start = 55;           % start frequency (Hz)
f_stop = 22000;         % stop frequency (Hz)
block_size = 2^12;      % block size (FFT window length)
ampl=0.1;                 % select peak amplitude of output re full scale

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


%% Generate a multi-sine output signal
% The number of output frequencies might be less than N_freqs, as the function removes redundant frequencies  
[t, sig, ms_indices] = generate_multisine(N_freqs, f_start, f_stop, fs, block_size, T);
sig=sig*ampl;  % scale signal
% ramp-up, ramp-down with hanning window
t_ramp = 200e-3;
n_ramp = floor(t_ramp*fs);
hann = hanning(2*n_ramp, 'periodic');
sig(1:n_ramp) = sig(1:n_ramp) .* hann(1:n_ramp);
sig(end-n_ramp+1:end) = sig(end-n_ramp+1:end) .* hann(n_ramp+1:end);

%% Play and record

[playData, recData, N_underrun, N_overrun] = play_rec(aPR, sig);
release(aPR);

t_recData = (0:size(recData, 1)-1)'/fs;

%% Process measurement data 
n_start = floor(t_start*fs);
n_stop = floor(t_stop*fs);
num_avg = floor((n_stop-n_start+1) /block_size); % number of blocks in data
n_stop = n_start + num_avg*block_size;          % align end to block border

rec = recData(n_start:n_stop-1,:);  % cut out data for analysis
t_rec = (0:length(rec)-1)/fs;

fprintf('Analyse data from T=%.2f to %.2f (averaging over %i blocks)\n',t_start,t_stop,num_avg);

% Average in time domain
rec_avg=mean(reshape(rec,block_size,num_avg,2),2);

%% FFT of input and output
% fft of the input
fft_rec = fft(rec_avg)/block_size;
fft_rec = 2*fft_rec(2:floor(block_size/2)+1,:); % single sided spectrum without DC
fft_freqs = fs*(1:block_size/2)'/block_size; % frequencies of spectrum

%% Spectrum for frequencies, where energy was provided
ms_fft_freqs = fft_freqs(ms_indices,:); % select frequency vector
ms_fft_rec = fft_rec(ms_indices,:);     % select frequencies in spectrum

% Save the averaged spectra
save('averaged_spectra.mat', 'ms_fft_rec', 'ms_fft_freqs');

% Load the data from the saved file
load('averaged_spectra.mat');

% Given that you have a known resistor value Rtest = 1 kΩ
Rtest = 1e3; % 1 kΩ = 1000 Ω

% U1 and U2 are the magnitudes of ms_fft_rec for each channel
U1 = abs(ms_fft_rec(:, 1));
U2 = abs(ms_fft_rec(:, 2));

% Calculate the impedance
Zin = U2 ./ (U1 - U2) * Rtest;

% Calculate the mean impedance
Zin_mean = mean(Zin);

% Plot the impedance over frequency
figure;
semilogx(ms_fft_freqs, abs(Zin));
title('Impedance over Frequency', 'FontSize',16);
xlabel('Frequency (Hz)', 'FontSize', 14);
ylabel('Impedance (Ohms)', 'FontSize',14);

% Display the mean impedance
fprintf('Mean Impedance: %f Ohms\n', Zin_mean);
