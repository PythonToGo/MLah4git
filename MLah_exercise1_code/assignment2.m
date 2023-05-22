clear all; clc; close all;

[y, fs] = audioread("tone_in_noise.wav");
t = linspace(0,length(y)/fs, length(y));
L = length(y);

% determine sample number 
tnew = [0.7 1.7]

samples = linspace(0.7, 1.7, length(y));

% cut off
figure(1);
plot(t,y);
xlim([0.7 1.7])





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


%{
title('Plotting sound file','FontSize',12);
xlabel('Time in s','FontSize',10,'FontWeight','normal');
ylabel('Amplitutde in Pa','FontSize',10,'FontWeight','normal');
%}