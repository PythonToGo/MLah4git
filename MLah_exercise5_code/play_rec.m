%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab-Subroutine for Measurement Laboratory at Home
% play_rec: writes and reads data from audio interface using aPR
% aPR is the MATLAB audioPlayerRecorder
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
% NOTE: the audio interface need an ASIO4all driver!!!
% NOTE: your set-up might need different added samples (line 14 and 16)!!
% by Ali Saeedi & Werner Hemmert, TUM, 20. Sept. 2020 - 27. Nov. 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [playData, recData, N_underrun,N_overrun] = play_rec(aPR, input)

buffer_size=aPR.BufferSize;        % read buffersize from aPR
added_samples_begin = 0* buffer_size; % pad extra block before signal start
% Replace 0 with 1 in some settings to prevent a strange blank area in the input. Test 0 and 1 to see this effect, at least in some configurations.
added_samples_end =4*buffer_size;     % pad extra blocks after the signal because of delay, make sure, you add enough buffers to read the full signal!

% prepare the input for two channels
if size(input, 2) == 1
    input = [input input];
end
N=size(input,1);
playData = [zeros(added_samples_begin, 2); ...
      input; zeros(ceil(N/buffer_size)*buffer_size-N, 2);...
      zeros(added_samples_end, 2)];  % construct output signal with multiple of block_size samples and add zeros before and after


recData = zeros(size(playData, 1), 2);  % recordings are saved in recData, allocate memory first!

n_total_buffers = floor(size(playData, 1)/buffer_size);

N_underrun=zeros(n_total_buffers,1);   % indicates buffer underrun 
N_overrun=N_underrun;                  % indicates buffer overrun
fprintf('recording...')
for i = 1:n_total_buffers              % iterate over number of buffers
    [recData(1+(i-1)*buffer_size:i*buffer_size, :),N_underrun(i),N_overrun(i) ]= aPR(playData(1+(i-1)*buffer_size:i*buffer_size, :));
end
fprintf('done\n')
%   Thats it already !
