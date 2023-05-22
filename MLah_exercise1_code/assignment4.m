clear; clc; close all;

% setting the variables
fs = 10240; 
T = 1/fs;
N = 1024;   
t = (0:N-1)/fs;
u1 = 0.5 + 0.1*sin(2*pi*100*t);
u2 = 0.1*sin(2*pi*105*t);

Y = fft(u2);
P2 = abs(Y/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(N/2))/N;



figure(401)
plot(t,u2,'r-');
title('u2(t) and its rms value','FontSize',12)
xlabel('time(t)');
ylabel('u2(t)');

hold on;
y = rms(u2,"all")       % y = 0.0707
plot(y,'bo');
hold off;


figure(402)
subplot(2,1,1);
plot(f, P1,'r-');
title('amplitude spectrum of u2(t)','FontSize',12)
xlabel('f (Hz)');
ylabel('|P1(f)|');

subplot(2,1,2)
loglog(f,P1,'g-');
title('u2(t) on log scale','FontSize',12)
xlabel('u2(t) on log');
ylabel('base-10 log scale');



%}
