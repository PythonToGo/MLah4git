clc
clear all
close all

% setting the variables
fs = 10240; 
T = 1/fs;
N = 1024;   
t = (0:N-1)/fs;
u2 = 0.1*sin(2*pi*105*t);


w_f = flattopwin(N)




% hann window
L = N+1;
w_h = hann(L);
wvtool(hann(L));
