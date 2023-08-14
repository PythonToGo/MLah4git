clc, close all, clear all;

T = 5;            % Total signal duration in (s), more averages with longer signal
t_start = 0.5;      % signal start: remove also transient part
t_stop=4.8;         % signal end: remove the fading part


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
%t_avg = (0:block_size-1)'/fs;  % Time vector for the averaged signal


%% FFT 
%%%%%%%%%%%%%%% A(b)%%%%%%%%%%%%%%%%
% U1 FFT
fft_U1 = fft(rec(:, 1)); % fft do
fft_U1 = fft_U1(1:block_size/2); % first half
frequencies = fs*(0:(block_size/2-1))/block_size; % frequency

% FFT of U2
fft_U2 = fft(rec(:, 2)); % fft do
fft_U2 = fft_U2(1:block_size/2); % first half

%% Plots
figure(1);
% U1
subplot(2,1,1)
plot(frequencies, 20*log10(abs(fft_U1)));
title('Amplitude spectrum U1(f)','FontSize',16);
xlabel('Frequency [Hz]','FontSize',14);
ylabel('Magnitude [dB]', 'FontSize',14);

subplot(2,1,2)
plot(frequencies, angle(fft_U1)*180/pi);
title('Phase spectrum U1(f)','FontSize',16);
xlabel('Frequency [Hz]','FontSize',14);
ylabel('Phase [degrees]','FontSize',14);


figure(2);
% U2
subplot(2,1,1)
plot(frequencies, 20*log10(abs(fft_U2)));
title('Amplitude spectrum U2(f)','FontSize',16);
xlabel('Frequency [Hz]','FontSize',14);
ylabel('Magnitude [dB]','FontSize',14);

subplot(2,1,2)
plot(frequencies, angle(fft_U2)*180/pi);
title('Phase spectrum U2(f)','FontSize',16);
xlabel('Frequency [Hz]','FontSize',14);
ylabel('Phase [degrees]','FontSize',14);


%% Calculate currents and voltages
%%%%%%%%% A(c) %%%%%%%%%%%%%%%%
Rm = 33000; % test resistancce

% current Iz to finden through Zt+Rm
Iz = fft_U1 / Rm; 

% voltage Uz
Uz = fft_U2 - fft_U1;

% impedance Zt
Zt = Uz ./ Iz;

% Plot impedance
% Magnitude
figure(3);
subplot(2,1,1)
plot(frequencies, 20*log10(abs(Zt)));
title('Impedance magnitude Zt(f)','FontSize',16);
xlabel('Frequency [Hz]','FontSize',14);
ylabel('Magnitude [dB]','FontSize',14);

% Phase
subplot(2,1,2)
plot(frequencies, angle(Zt)*180/pi);
title('Impedance phase Zt(f)','FontSize',16);
xlabel('Frequency [Hz]','FontSize',14);
ylabel('Phase [degrees]','FontSize',14);

%% Plot real and imaginary parts of Zt
%%%%%%%%% 1(d) %%%%%%%

figure(4);
% Real part
subplot(2,1,1);
plot(frequencies, real(Zt));
title('Real part of Zt(f)','FontSize',16);
xlabel('Frequency [Hz]','FontSize',14);
ylabel('Real [Ohm]','FontSize',14);

% Imaginary part
subplot(2,1,2);
plot(frequencies, imag(Zt));
title('Imaginary part of Zt(f)','FontSize',16);
xlabel('Frequency [Hz]','FontSize',14);
ylabel('Imaginary [Ohm]','FontSize',14);


%% Calculate the total admittance
%%%%%%%%%%%%%% 1(e) %%%%%%%%%%%%
Yt = 1 ./ Zt;

% Conductance Gt, Susceptance Bt = Real(Yt), Imag(Yt)
Gt = real(Yt);
Bt = imag(Yt);

% Rt, Ct
Rt = 1 ./ Gt;
Ct = Bt ./ (2*pi*frequencies);

% Now, let's plot Rt and Ct

% plot Rt, Ct
figure(5);
% Rt
subplot(2,1,1);
semilogx(frequencies, Rt);
title('Resistor Rt(f)','FontSize',16);
xlabel('Frequency [Hz]','FontSize',14);
ylabel('Resistance [Ohm]','FontSize',14);

% Ct
subplot(2,1,2);
semilogx(frequencies, Ct);
title('Capacitor Ct(f)','FontSize',16);
xlabel('Frequency [Hz]','FontSize',14);
ylabel('Capacitance [F]','FontSize',14);
 
%% get average over frequency
%%%%%%%%%% 1(f) %%%%%%%%%%%%
Rt_avg = mean(Rt);
Ct_avg = mean(Ct);

fprintf('The average value of Rt over all frequencies is %.2f Ohm.\n', Rt_avg);
fprintf('The average value of Ct over all frequencies is %.2e F.\n', Ct_avg);

