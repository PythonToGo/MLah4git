clear; clc; close all;

% setting the variables
fs = 10240; % sampling frequency
T = 1/fs;
N = 1024;


t = (0:N)/fs;
u1 = 0.5 + 0.1*sin(2*pi*100*t);
Y = fft(u1);
P2 = abs(Y/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(N/2))/N;

% tnew
tnew = (0:N-1)/fs;
u1new = 0.5 + 0.1*sin(2*pi*100*tnew);
Ynew = fft(u1new);
P2new = abs(Ynew/N);
P1new = P2new(1:N/2+1);
P1new(2:end-1) = 2*P1new(2:end-1);
fnew = fs*(0:(N/2))/N;


figure(301)

%u1(t) with t=(0:N)/fs
subplot(3,2,1)
plot(t, u1,'r-');
title('u1(t)','FontSize',12);
xlabel('time(t)');
ylabel('u1(t)');

% u1(tnew) 
subplot(3,2,2)
plot(tnew, u1new,'b-');
title('u1(tnew)','FontSize',12)
xlabel('time(tnew)');
ylabel('u1(tnew)');


% SSA of u1(t) 
subplot(3,2,3);
plot(f,P1,'r-');
title('Single-Sided Amp spec of u1(t)','FontSize',12)
xlabel('f (Hz)');
ylabel('|P1(f)|');

% SSA of u1(tnew) 
subplot(3,2,4)
plot(fnew,P1new,'b-');
title('Single-Sided Amp Spec of u1(tnew)','FontSize',12)
xlabel('fnew (Hz)');
ylabel('|P1(fnew)|');


subplot(3,2,5)
loglog(f,P1,'r-');
title('u1(t) on a log-log scale','FontSize',12)
xlabel('u1(t) on log');
ylabel('base-10 log scale');

subplot(3,2,6);
loglog(fnew,P1new,'b-');
title('u1(tnew) on a log-log scale','FontSize',12)
xlabel('u1(tnew) on log');
ylabel('base-10 log scale');


sgtitle('    Left side t=(0:N)/fs                         Right side tnew=(0:N-1)/fs')


