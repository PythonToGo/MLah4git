clear; clc; close all;

% setting the variables

fs = 48000; 
T= 0.5;
N = fs*T;
t = (0:N-1)/fs;
u3 = 0.1*cos(2*pi*230*t);

plot(t,u3);
title('Plotting u3(t)','FontSize',12);
xlabel('Time in s','FontSize',10,'FontWeight','normal');
ylabel('u(t)','FontSize',10,'FontWeight','normal');

sound(u3,fs);

