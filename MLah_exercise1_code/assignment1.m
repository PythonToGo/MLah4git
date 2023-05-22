clear all; clc; close all;

%file read
[y, fs] = audioread("tone_in_noise.wav");

% time domain
t = linspace(0,length(y)/fs, length(y));
L = length(y);


% assignment 1a)
figure(1);
plot(t,y,'b-');
title('Plotting sound file','FontSize',12);
xlabel('Time in s','FontSize',10,'FontWeight','normal');
ylabel('Amplitutde in Pa','FontSize',10,'FontWeight','normal');

% FFT change
Y = fft(y);
P2 = abs(Y/L);  % double-sided
P1 = P2(1:L/2+1);   % single-sided
P1(2:end-1) = 2*P1(2:end-1);


%for plotting
figure(2);
freq = fs*(0:(L/2))/L;
plot(freq, P1, 'ro');
title('Singsle-Sided Amplitude Spectrum')
xlabel('Frequency(Hz)','FontSize',12,'FontWeight','normal')
ylabel('Magnitude','FontSize',10,'FontWeight','normal')


% using NFFT
NFFT = 2^nextpow2(L);
Y = fft(y, NFFT)/L;
f = fs/2*linspace(0,1,NFFT/2);


figure(3);
plot(f,2*abs(Y(1:NFFT/2)));
title('Singsle-Sided Amplitude Spectrum using NFFT')
xlabel('Frequency(Hz)','FontSize',12,'FontWeight','normal')
ylabel('Magnitude','FontSize',10,'FontWeight','normal')


% plotting it using loglog
figure(4);
loglog(f,2*abs(Y(1:NFFT/2)),'r-');
title('plotting single-sided amplitude spectrum using loglog')
xlabel('Frequency(Hz)','FontSize',12,'FontWeight','normal')
ylabel('Magnitude','FontSize',10,'FontWeight','normal')



