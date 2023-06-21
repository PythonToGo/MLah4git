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

% figure;
% plot(t, sig);
% xlabel('Time in s');
% label('Amplitude re full scale');

%% Play and record

[playData, recData, N_underrun, N_overrun] = play_rec(aPR, sig);
release(aPR);

t_recData = (0:size(recData, 1)-1)'/fs;

%% Plot raw data
figure;
subplot(211); 
plot(t_recData, playData);
sgtitle('Played & recorded data');
title('Played data')
legend('Ch 1', 'Ch 2');
xlabel('Time (s)');
ylabel('Amplitude re full scale');

subplot(212); 
plot(t_recData, recData);
title('Recorded data (raw)');
legend('Ch 1', 'Ch 2');
xlabel('Time (s)');
ylabel('Amplitude re full scale');

%% Process measurement data 
n_start = floor(t_start*fs);
n_stop = floor(t_stop*fs);
num_avg = floor((n_stop-n_start+1) /block_size); % number of blocks in data
n_stop = n_start + num_avg*block_size;          % align end to block border

rec = recData(n_start:n_stop-1,:);  % cut out data for analysis
t_rec = (0:length(rec)-1)/fs;

fprintf('Analyse data from T=%.2f to %.2f (averaging over %i blocks)\n',t_start,t_stop,num_avg);
%% Plot data selected for analysis into raw data
figure
subplot(2,1,1)
sgtitle('Data selected for analysis');
plot(t_recData, recData(:,1),'k');
hold on
plot(t_rec+n_start/fs, rec(:, 1), 'b');
ylabel('CH1 Amplitude re full scale');
subplot(2,1,2)
plot(t_recData, recData(:,2),'k');
hold on
plot(t_rec+n_start/fs, rec(:, 2), 'r');
xlabel('Time/s');
ylabel('CH2 Amplitude re full scale');

% Average in time domain
rec_avg=mean(reshape(rec,block_size,num_avg,2),2);

%% FFT of input and output
% fft of the input
fft_rec = fft(rec_avg)/block_size;
fft_rec = 2*fft_rec(2:floor(block_size/2)+1,:); % single sided spectrum without DC
fft_freqs = fs*(1:block_size/2)'/block_size; % frequencies of spectrum
    

%% plot raw spectrum
figure
subplot(2,1,1)
sgtitle('Spectrum of data selected for analysis');
stem(fft_freqs/1000, abs(fft_rec(:, 1)), 'k');
set(gca,'XScale','log','YScale','log');
ylabel('CH1 Magnitude re full scale');
subplot(2,1,2)
stem(fft_freqs/1000, abs(fft_rec(:, 2)), 'r');
set(gca,'XScale','log','YScale','log');
xlabel('f/kHz');
ylabel('CH2 Magnitude re full scale');

%% Spectrum for frequencies, where energy was provided
ms_fft_freqs = fft_freqs(ms_indices,:); % select frequency vector
ms_fft_rec = fft_rec(ms_indices,:);     % select frequencies in spectrum

figure
subplot(2,1,1)
loglog(ms_fft_freqs, abs(ms_fft_rec(:, 1)), 'k.');
ylabel('CH1 Magnitude re full scale');
%axis([10 fs/2 1e-8 1])
subplot(2,1,2)
loglog(ms_fft_freqs, abs(ms_fft_rec(:, 2)), 'r.');
xlabel('f/Hz');
ylabel('CH2 Magnitude re full scale');
sgtitle('CH1 & CH2 Spectra')
% axis([10 fs/2 1e-8 1])

%%
% numerator and denominator
num = [100000000];
den = [1 30000 100000000];

% transfer function
sys = tf(num, den);

% Generate the Bode plot
bode(sys);
grid on;

%% Calculate transfer function
H=ms_fft_rec(:, 2) ./ ms_fft_rec(:, 1);
figure
subplot(2,2,1)
sgtitle('Compare calculated Transfer Function with theoretical');
loglog(ms_fft_freqs, abs(H), 'k.');
ylabel('|H|');
subplot(2,2,3)
semilogx(ms_fft_freqs, unwrap(angle(H))*180/pi, 'k.');
xlabel('f/Hz');
ylabel('Phase(H)/deg');

subplot(2,2,[2,4])
bode(sys);
grid on;

%% --------------------------------- plot --------------------------------
figure
fontSize=8;
subplot(2,1,1);        % create plots for magnitude and phase 
loglog(ms_fft_freqs/1000, abs(H), 'k.');
%xlabel('Frequenz / kHz','FontSize',fontSize) %CAREFUL: Plot scaled in kHz
ylabel('|H(f)|','FontSize',fontSize)
axis([50/1000 24 0.5 2])                       % scale plot

% ------------- some tricks for plotting --------------------------------
%y_pos=[0.05 0.1 0.2 0.5 1];                   % position for y-axes labels
%  set(gca,'YtickLabel',y_pos,'FontSize',fontSize);
%set(gca,'YTick',y_pos)
x_pos=[0.1  1  10 100];      % define labels for x-axes
set(gca,'XtickLabel',x_pos,'FontSize',fontSize);
% --------------------------- plot phase --------------------------------
subplot(2,1,2);                               % plot phase below
semilogx(ms_fft_freqs/1000, unwrap(angle(H),pi/2)*180/pi, 'k.');
semilogx(ms_fft_freqs/1000, angle(H)*180/pi, 'k.');
xlabel('f in kHz','FontSize',fontSize);
ylabel('Phase(H) in deg','FontSize',fontSize);
axis([50/1000 24 -370 100])                       % scale plot
y_pos=[-360 -270 -180 -90 0 90];                  % position for y-axes labels
set(gca,'YTick',y_pos)
set(gca,'YtickLabel',y_pos,'FontSize',fontSize);
x_pos=[0.1  1  10 100];      % define labels for x-axes
set(gca,'XtickLabel',x_pos,'FontSize',fontSize);
%  H=line([0.01 0.1],[0 0]);                  % Bode limits
%  set(H,'LineStyle','-.','Color','r')
%  H=line([0.1 10],[0 -90]);                  % Bode limits
%  set(H,'LineStyle','-.','Color','r')
%  H=line([10 20],[-90 -90]);                 % Bode limits
%  set(H,'LineStyle','-.','Color','r')

subplot(2,1,1);        % create plots for magnitude and phase 
title(name); name;
sgtitle('Transfer Function Plot')
% set(gcf,'Units','Centimeters','Position',[0 1 8.4 9],'PaperPositionMode','auto')
print(name, '-depsc')                      % create scaleable figure
print(name, '-dtiff', '-r300')             % cretes pixel figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------- D O N E -----------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
