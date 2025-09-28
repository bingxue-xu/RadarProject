%% Range-Doppler Measurements - FMCW Mode
% 
% Group 7
% Review Meeting 3
% 23/09/2025
% 

%% Cleaning
clear; 
clc;

%% Parameters
c = 3e8; % speed of light (m/s)
f_start = 2.405e9; % start freq (Hz)
f_stop = 2.495e9; % stop freq (Hz)
deltaf = f_stop - f_start; % sweep bandwidth (Hz)
Tp = 5e-3; % chirp duration (s) 
fc = 0.5*(f_start + f_stop); % center frequency (Hz)
lambda = c / fc; % wavelength (m)

%% Read data
%our data
[data, Fs] = audioread("fmcw_18_t4_2targets.wav");
radar_backscatter = data(:,2);
sync_signal = data(:,1);

%prof data
%[data, Fs] = audioread("Range_Test_File2.m4a");
%radar_backscatter = data(:,1);
%sync_signal = data(:,2);

shift_amount = floor(Tp*Fs*0.5);  % number of indices to shift
sync_signal = [sync_signal(shift_amount+1:end); zeros(shift_amount,1)];


chunktime = 1; % s For subtask: changing the length of the chunk of data being processed
chunk = Fs * chunktime;
Indx = 1;
threshold = 0;

%% Video Output
v = VideoWriter('Range_Doppler_Spectrogram_MS_MTI.avi');
v.FrameRate = 5;
open(v);
fig = figure('Color','w');


%% Chunk processing loop
while Indx + chunk - 1 <= length(radar_backscatter)

    rb_chunk = radar_backscatter(Indx : Indx + chunk - 1);
    sync_chunk = sync_signal(Indx : Indx + chunk - 1);

    syncPulse = zeros(size(sync_chunk));
    syncPulse(sync_chunk > threshold) = 1;

    d_sync = diff([0; syncPulse; 0]); % add zeros incase mid up/down chirp for detection
    start_idx = find(d_sync == 1);
    end_idx = find(d_sync == -1) - 1;

    N_chirp = round(Fs * Tp); % fixed chirp length
    num_chirps = numel(start_idx);

    upchirp_matrix = zeros(num_chirps, N_chirp); % initializing the matrix with 0s
    
    % using METHOD 2 (reference to SG2_range_fmcw.m): if the up chirp is not complete, remain data is filled with 0
    for k = 1:num_chirps
        if end_idx(k) - start_idx(k) + 1 >= N_chirp
            upchirp_matrix(k,:) = rb_chunk(start_idx(k) : start_idx(k) + N_chirp - 1);
        else
            % if a chirp is shorter than expected -> 0-pad
            seg = rb_chunk(start_idx(k):end_idx(k));
            upchirp_matrix(k,1:numel(seg)) = seg;
        end
    end

    % Raw data - no CR
    data_raw = upchirp_matrix;

    % MS clutter rejection
    upchirp_ms = upchirp_matrix - mean(upchirp_matrix, 1);

    % 2-pulse MTI
    data_2MTI = upchirp_ms(2:end, :) - upchirp_ms(1:end-1, :);

    %IFFT
    Y1raw = ifft(data_raw, [], 2); % fast-time IFFT for range
    Y2raw = ifft(Y1raw, [], 1); % slow-time IFFT for velocity
    Y2absraw = abs(Y2raw);
    Y2_dBraw = 20*log10(Y2absraw);

    Y1ms = ifft(upchirp_ms, [], 2); % fast-time IFFT for range
    Y2ms = ifft(Y1ms, [], 1); % slow-time IFFT for velocity
    Y2absms = abs(Y2ms);
    Y2_dBms = 20*log10(Y2absms);

    N = length(upchirp_ms);

    %this is part thats used i think
    Y1 = ifft(upchirp_ms, [], 2); % fast-time IFFT for range %CHANGED: use ms no MTI
    Y2 = ifft(Y1, [], 1); % slow-time IFFT for velocity
    %Y2 = ifftshift(Y2, 1); %didnt do anything useful
    Y2abs = abs(Y2);
    Y2_dB = 20*log10(Y2abs);

    Y1_zp = ifft(data_2MTI, 4*N, 2); % fast-time IFFT for range
    Y2_zp = ifft(Y1_zp, 4*N, 1); % slow-time IFFT for velocity
    Y2abs_zp = abs(Y2_zp);
    Y2_dB_zp = 20*log10(Y2abs_zp);


    Y2_lower_raw = Y2_dBraw(:, 1:floor(N/2));
    Y2_lower_raw = Y2_lower_raw - max(Y2_lower_raw(:));

    Y2_lower_ms = Y2_dBms(:, 1:floor(N/2));
    Y2_lower_ms = Y2_lower_ms - max(Y2_lower_ms(:));

    Y2_lower = Y2_dB(:, 1:floor(N/2)); %
    Y2_lower = Y2_lower - max(Y2_lower(:));

    Y2_lower_zp = Y2_dB_zp(:, 1:floor(4*N/2));
    Y2_lower_zp = Y2_lower_zp - max(Y2_lower_zp(:));

    %% No Zero Padding
    % Velocity axis
    N_slow = size(data_2MTI, 1);
    df_doppler = 1 / (Tp * N_slow);            
    velocity = (lambda/2) * (0:(floor(N_slow/2)-1)) * df_doppler; % positive Doppler

    % Range axis using deltaR and Rmax
    deltaR = c / (2 * deltaf);
    Rmax = (N * deltaR) / 2;
    range = linspace(0, Rmax, floor(N/2));

    %% Zero Padding
    % B = f_stop - f_start;
    % slope = B / Tp; % Hz/s - chirp slope
    % 
    % % Range axis (from fast-time FFT bins)
    % df_range = Fs / (4*N);
    % f_range = (0:(4*N/2 - 1)) * df_range;
    % range_axis = (c * f_range) / (2 * slope);
    % range = range_axis; % meters
    % 
    % % Velocity axis (from slow-time FFT bins)
    % N_doppler = size(Y2_zp,1);
    % df_dopp = (1 / Tp) / N_doppler;
    % f_dopp = (1 : N_doppler/2-1) * df_dopp;
    % velocity = (lambda/2) * f_dopp; % m/s

    %% Plot frame
    imagesc(range, velocity, Y2_lower, [-45 0]); 
    % imagesc(range, velocity, Y2_lower_zp, [-50 0]); % for zero padding
    axis ij;
    xlabel('Range (m)'); ylabel('Velocity (m/s)');
    title('Range–Doppler Plot (IFFT–IFFT, with MS and 2-pulse MTI)');
    % title('Range–Doppler Plot (IFFT–IFFT, with MS and 2-pulse MTI) zero padding (4N)'); % for zero padding
    colorbar;
    xlim([0 100]);
    ylim([0 6]);

    drawnow;
    frame = getframe(fig);
    writeVideo(v, frame);

    % Increment the index by the chunk amount
    Indx = Indx + chunk;
end

close(v);



