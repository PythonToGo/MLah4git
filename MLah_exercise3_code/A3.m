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
t_stop=8;           % signal analysis end: remove the fading part

f = 100;            % set signal frequency in Hz
N = floor(T*fs);    % calculate number of samples
t = (0:N-1)/fs;    % generate time vector in s
amp = 1.0;          % set peak amplitude of output re full scale <=1


%  white noise
white_noise = rand(N, 1) * 0.5 - 0.25;   
white_noise(:,2) = rand(N,1) * 0.5 - 0.25;    

% ramp-up, ramp-down signal with hanning window
t_ramp = 0.1;                    % ramp-up and down time in s
n_ramp = floor(t_ramp*fs);          % calculate number of samples for ramp
hann = hanning(2*n_ramp, 'periodic');     % ramp is a hanning window shape
hann(:,2) = hanning(2*n_ramp, 'periodic');
white_noise(1:n_ramp,:) = white_noise(1:n_ramp,:) .* hann(1:n_ramp,:); % apply ramp on
white_noise(end-n_ramp+1:end,:) = white_noise(end-n_ramp+1:end,:) .* hann(n_ramp+1:end,:);

%
%sig = amp*sin(2*pi*f*t);            % calculate sine signal
%sig(:,2) = amp*sin(2*pi*f*t);       % set signal for second channel

% Play and Record signal
fprintf('Start recording  -- read output voltage with voltmeter (RMS)\n...\n')
[playData, recData] = play_rec(aPR, white_noise);
fprintf('End recording\n')
t_recData = (0:length(recData)-1)/fs;

release(aPR);
%rms_amp=input('RMS amplitude in V_rms:');    % input RMS value
rms_amp = 0.5;
fprintf('this coresponds to a peak amplitude of %.4f V\n',rms_amp*sqrt(2));

figure(10);
plot(t_recData, recData);
xlabel('Time (s)','FontSize',14);
ylabel('Amplitude','FontSize',14);
title('Original recorded Signal','FontSize',16);

%% Assignment 3(a)
u1t = recData((145000:147000), 1);  % Save input signal u1(t) to u1t variable
u2t = recData((145000:147000), 2);  % Save input signal u2(t) to u2t variable
t_recData_seg = t_recData(145000:147000);

% Plot Signal u1(t) and u2(t)
figure(1);
% u1(t)
subplot(2, 1, 1);
plot(t_recData_seg, u1t);
xlabel('Time (s)','FontSize',14);
ylabel('Amplitude','FontSize',14);
title('Input Signal u1(t)','FontSize',16);
% u2(t)
subplot(2, 1, 2);
plot(t_recData_seg, u2t);
xlabel('Time (s)','FontSize',14);
ylabel('Amplitude','FontSize',14);
title('Output Signal u2(t)','FontSize',16);


%% Assignment 3(b)

% length, freqency
L = length(u1t);  % Length of the signals
f_sss = (0:L/2)*(fs/L);  % Frequency vector for single-sided spectrum

% Fourier Transform of the u1(t) / u2(t)
U1f = fft(u1t);
U2f = fft(u2t);


% single sided amplitude spectrum of U1f / U2f
U1f_single = abs(U1f(1:L/2+1))/L;
U2f_single = abs(U2f(1:L/2+1))/L;


%% Assingment 3(b)
% Plot single-sided amplitude spectra
figure(2);
subplot(2, 1, 1);
plot(f_sss, U1f_single);
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Amplitude','FontSize',14);
title('Single-Sided Amplitude Spectrum of U1(f)','FontSize',16);

subplot(2, 1, 2);
plot(f_sss, U1f_single);
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Amplitude','FontSize',14);
title('Single-Sided Amplitude Spectrum of U2(f)','FontSize',16);



%% 
% Transfer function H(f)
Hf = U2f ./ U1f;


figure(1);
plot(Hf);
title('Transfer Function H(f)', 'FontSize',16);


% Frequency vector
%f_bode = logspace(0, log10(fs/2), length(Hf));
%[mag,phase] = bode(Hf,f_bode);

% Frequency = vector form; 0Hz ~24kHz
%-------------------------------QEUSTION for # fs vector_________________
f_bode = linspace(0, fs, length(Hf)); 


mag = abs(Hf);


%%
% magnitude
mag = abs(Hf);
% Create the Bode plot
figure;
subplot(2,1,1);
semilogx(f_bode, 20*log10(mag));
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Magnitude (dB)','FontSize',14);
title('Magnitude Response','FontSize',16);
grid on;

subplot(2,1,2);
semilogx(f_bode, unwrap(angle(Hf))*180/pi);
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Phase (degrees)','FontSize',14);
title('Phase Response','FontSize',16);
grid on;

%%
