% Define the signal u2(t)
fs = 1000; % Sampling frequency
t = (0:1/fs:1-1/fs); % Time vector
u2 = 0.1*sin(2*pi*105*t); % Signal u2(t)

% Apply a flat top window to the signal
win = flattopwin(length(u2))'; % Flat top window
u2_win = u2 .* win; % Signal with flat top window applied

% Calculate the single-sided amplitude spectrum of the signal
N = length(u2_win);
U2 = fft(u2_win)/N;
U2_amp = abs(U2(1:N/2+1));
U2_amp(2:end-1) = 2*U2_amp(2:end-1);
N = length(u2); % length of the signal
Y = fft(u2) / N; % calculate the Fourier coefficients and normalize
f = (0:N-1) * fs / N; % calculate the frequency vector
spectrum = 2 * abs(Y(1:N/2)); % extract the single-sided spectrum


% Calculate the window correction factor and CPG
A_window = sum(win)/length(win);
correction_factor_ft = sqrt(1/A_window)
cf = 1/A_window;
CPG = 20*log10(cf*U2_amp);


% Plot the spectrum on a log-log scale
f = (0:N/2)*fs/N;
loglog(f, U2_amp);
xlabel('Frequency (Hz)');
ylabel('Amplitude (V)');
title('Single-sided Amplitude Spectrum of u2(t) with Flat Top Window');


% Calculate the amplitude of the tone
max_spectrum = max(abs(spectrum)); % find maximum amplitude in spectrum
amplitude = max_spectrum / correction_factor_ft; % divide by correction factor

% Display the amplitude and CPG
disp(['The amplitude with a float top window is ', num2str(amplitude)])
disp(['The CPG with a flat top window is ', num2str(CPG)])

% Plot the CPG
figure(30)
plot(CPG);
title('CPG for a flat top window of u2(t)');


% Define the Hanning window
win = hann(length(u2));

% Calculate the correction factor with hanning window
correction_factor_hann = sum(win.^2)/length(win);
disp(['The correction factor with' ...
    'hanning window is ', num2str(correction_factor_hann)])

% calculate CPG in dB with hannind window
CPG_hann = mag2db(max(U2_amp)/rms(U2_amp)); 
disp(['The CPG with hanning window is ',...
    num2str(CPG_hann)])
