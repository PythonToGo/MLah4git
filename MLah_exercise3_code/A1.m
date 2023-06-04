clear; clc; close all;

%% RC CIRCUIT SETTING 
R = 2000; % 오오오오옴
C = 100e-9; % farad unit changed.

% set frequency f :  0 ~ 24kHz
f = logspace(0, log10(24000), 1000);

% set s and omega :  s = jw , w
w = 2*pi*f; % Converting frequency -> angular frequency
s = 2i*pi*f;  

%% Trnasfer function
H = @(s)(s)./(s+5000);

% magnitude estimation
mag = abs(H(1j*w));
magDB = 20 * log10(mag);
% phase estimation
phaseDeg = rad2deg( angle(H(1j*w)) );


%% Assignment 1(a) -> Report
%% Assignment 1(b) -> Report
%% Bode plot & Transfer Function plot :: Assignment 1(c)


% plotting
figure(2);
subplot(211);
  semilogx(w/2/pi,magDB);
  xlabel('frequency [Hz]','FontSize', 14);
  ylabel('magniude [dB]','FontSize', 14);
  title('Bode diagram','FontSize', 16);
  grid on;
subplot(212);
  semilogx(w/2/pi,phaseDeg);
  xlabel('frequency [Hz]','FontSize', 14);
  ylabel('phase angle [°]','FontSize', 14);
  title('Phase diagram','FontSize', 16);
  grid on;


