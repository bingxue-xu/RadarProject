%% Real-time Measurements - FMCW Mode
% 
% Group 7
% Review Meeting 3
% 23/09/2025
% 


%% FMCW Radar – Real-time / Pseudo real-time Range Measurements (Batch of Chirps)
clear; clc; close all;

%% USER SETTINGS
mode = "real";   % choose: "real" or "pseudo"
filename = "Test_files/fmcw_18_t3_2targets.wav"; % only used in pseudo mode

Tp = 5e-3;             % chirp duration (s)
f_start = 2.4e9;       % Hz
f_stop  = 2.5e9;       % Hz
c = 3e8;               % m/s
threshold = 0;         % sync threshold
batchSize = 10;        % number of up-chirps per batch

%% ACQUISITION SETUP
if mode == "real"
    Fs = 44100; % must match soundcard
    devices = getAudioDevices(audioDeviceReader);
    device_number = 3; % CHOOSE BEFORE RUN
    device_name = devices{device_number};

    deviceReader = audioDeviceReader( ...
        'Device', device_name, ...
        'SampleRate', Fs, ...
        'NumChannels', 2, ...
        'SamplesPerFrame', round(Tp*Fs*batchSize*2) ); % at least batchSize chirps
    disp('Real-time FMCW acquisition started');
else
    [data, Fs] = audioread(filename);
    radar_data = data(:,2); % radar channel
    sync_data  = data(:,1); % sync channel
    shift_amount = floor(Tp*Fs*0.5);
    sync_data = [sync_data(shift_amount+1:end); zeros(shift_amount,1)];
    disp('Pseudo real-time playback started');
end

%% VIDEO SETUP
saveVideo = true;
if saveVideo
    v = VideoWriter('FMCW_videos/FMCW_RangeTime.avi');
    v.FrameRate = 10;
    open(v);
end

%% FFT PARAMETERS
N = round(Tp*Fs);       % samples per chirp
fft_size = 4*N;         
B = f_stop - f_start;   
slope = B/Tp;           
df = Fs/fft_size;       

Rmax = c*Fs/(2*slope);  % max unambiguous range
range_axis = (0:fft_size/2-1) * (c/(2*slope)) * df;

%% PLOT INITIALIZATION
fig = figure('Color','w');
hImg = imagesc(range_axis, 0, zeros(1,length(range_axis)), [-50 0]);
axis xy;
xlabel('Range (m)');
ylabel('Time (s)');
set(gca,'YDir','reverse');
title('FMCW Range–Time Plot');
colorbar;
xlim([0 15]);

%% ESC key stop control
stopFlag = false;
set(fig, 'KeyPressFcn', @(src, event) assignin('base','stopFlag', strcmp(event.Key,'escape')));

%% PROCESSING LOOP
frameCount = 0;
t_axis = [];
range_map = [];

%% Initialize chirp counter
totalChirps = 0;

t_start = tic;

while ~stopFlag
    % --- ACQUIRE FRAME ---
    if mode == "real"
        frameStereo = deviceReader();
        if isempty(frameStereo), continue; end
        radar_frame_full = frameStereo(:,2);
        sync_frame_full  = frameStereo(:,1);

        % detect up-chirp indices in the frame
        sync_up = sync_frame_full > threshold;
        d_sync = diff([0; sync_up; 0]);
        start_idx_frame = find(d_sync == 1);

        % collect batchSize chirps
        batch_matrix = zeros(batchSize, N);
        chirp_count = 0;
        for k = 1:numel(start_idx_frame)
            idx_s = start_idx_frame(k);
            if idx_s + N - 1 > length(radar_frame_full), break; end
            chirp_count = chirp_count + 1;
            batch_matrix(chirp_count,:) = radar_frame_full(idx_s:idx_s+N-1);
            if chirp_count == batchSize, break; end
        end
        if chirp_count < batchSize, continue; end  % wait for full batch

        % update total chirps
        totalChirps = totalChirps + batchSize;

    else
        % --- PSEUDO MODE ---
        if frameCount == 0
            sync_up = sync_data > threshold;
            d_sync = diff([0; sync_up; 0]);
            start_idx = find(d_sync == 1);
            M = numel(start_idx);
        end

        % collect batch
        batch_matrix = zeros(batchSize, N);
        chirp_count = 0;
        for k = 1:batchSize
            frameCount = frameCount + 1;
            if frameCount > M, break; end
            idx_s = start_idx(frameCount);
            if idx_s + N - 1 > length(radar_data), break; end
            chirp_count = chirp_count + 1;
            batch_matrix(chirp_count,:) = radar_data(idx_s:idx_s+N-1);
        end
        if chirp_count < batchSize, break; end

        % update total chirps
        totalChirps = totalChirps + batchSize;
    end

    % --- PROCESS BATCH ---
    batch_matrix = batch_matrix - mean(batch_matrix, 1);
    Y = fft(batch_matrix, fft_size, 2);
    Y_mag = abs(Y(:,1:fft_size/2));
    Y_dB = 20*log10(Y_mag);
    Y_norm = Y_dB - max(Y_dB(:));

    % --- UPDATE DATA STRUCTURES ---
    % t_axis = [t_axis; (totalChirps-batchSize+1:totalChirps)'*Tp]; %#ok
    t_now = toc(t_start);
    t_batch = t_now - Tp*(batchSize-1:-1:0)';
    t_axis = [t_axis; t_batch]; %#ok


    range_map = [range_map; Y_norm]; %#ok

    % --- UPDATE PLOT ---
    set(hImg, 'XData', range_axis, 'YData', t_axis, 'CData', range_map);
    ylim([-Inf t_axis(end)]);
    drawnow;

    % --- SAVE VIDEO ---
    if saveVideo
        frameObj = getframe(fig);
        writeVideo(v, frameObj);
    end
end


if saveVideo
    close(v);
    disp('Video saved successfully.');
end
