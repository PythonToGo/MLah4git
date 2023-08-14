%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab-Subroutine for Measurement Laboratory at Home
% generate_multisine: generates a multi-sine of N_freqs 
% log-spaced frequencies between f_start and f_stop with random phases.
% Please note that the signal will always have a maximum of N_freqs but
% less, when the linear FFT resolution is coarser than the log spaced
% resolution (usually happens at low frequencies).
%% 7.May 2021 replace randi with rand 
%% also add argument to initialize rng to get same sequence if needed
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
% NOTE: Signal duration should be >= of the block_size
% by Ali Saeedi & Werner Hemmert, TUM, 20. Sept. 2020 - 06. Dec. 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, sig, ms_indices] = generate_multisine(N_freqs, f_start, f_stop, fs, block_size, T, rng_seed)

if(nargin==7)
   rng(rng_seed);    % initialize random number generator if requested
end

if T < block_size/fs
    fprintf('Error: signal length is shorter than block size\n');
end
if f_stop > fs/2
    fprintf('Error: f_stop > fs/2  setting f_stop to fs/2\n');
    f_stop=fs/2
end
        
    f1 = log10(f_start);
    f2 = log10(f_stop);
    log_spaced_freqs = logspace(f1, f2, N_freqs)';	% N_freqs log-spaced frequencies between f_start and f_stop
    linear_spaced_freqs = fs*(1:block_size/2)'/block_size;	% linear-spaced frequencies (= FFT frequencies)
    
    % Normalize log-spaced frequencies to FFT frequency resolution. This helps to find the closest FFT frequencies to log-spaced frequencies as well as to have precise frequency bins.
    % This also helps to have precise frequeny bins.
    norm_log_freqs = round(log_spaced_freqs/(fs/2)*block_size/2);    

    % remove repetitive indices.
    ms_indices = unique(norm_log_freqs);
    
    % The final number of frequencies might be less than N_freqs, but sitll enough number of frequencies will remain.
    new_N = length(ms_indices);
    
    ms_freqs = linear_spaced_freqs(ms_indices);	% Frequencies of the multi-sine signal, omit DC!!

    amps = ones(new_N, 1);	% same normalized amplitude for all frequencies
    
    L_sig = floor(T*fs);
    t = (0:L_sig-1)'/fs;
    
    sig = zeros(L_sig, 1);
    % Generate random phases and convert them from degrees into radian.
    % Random phases prevent peak alignments
%%    rand_phases = (pi/180) * randi(360, new_N, 1);   %% Do NOT use INTEGERS!!
    rand_phases = 2*pi * rand(new_N, 1);
    
    for iF = 1:new_N
        
        % Generate single sine of each frequency
        temp_sine = amps(iF)*sin(2*pi*ms_freqs(iF)*t + rand_phases(iF));
        
        % Add them together
        sig = sig + temp_sine;        
        
    end
    
    % Find maximum of absolute value of the multi-sine signal
    max_sig = max([abs(min(sig)), abs(max(sig))]);
    
    % Normalize the signal amplitude to 1
    sig = sig/max_sig;
   
end

