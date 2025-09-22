%% Real Time Range Measurements - FMCW Radar
%% READ:
% 
% In pseudo real-time mode with file, wait for execution to finish so as 
% not to corrupt the video
% 
% In real-time mode:
% 1. Identify the correct input device to acquire the data
% 2. To finish acquisition/video, simply press the 'ESC' key to exit the 
% program correctly and avoid corrupting the video
% 
%% Cleaning
clear; clc; close all;

%% USER SETTINGS
mode = "pseudo";   % choose: "real" or "pseudo"
filename = "TEST_FILES18/fmcw_18_t1.wav"; % only used in pseudo mode

% constants
Tp = 5e-3;            % s
start_freq = 2.405e9; % Hz
stop_freq  = 2.495e9; % Hz
threshold = 0;        % sync threshold
c = 299792458;        % m/s

%% ACQUISITION SETUP
if mode == "real"
    % live from soundcard (Audio Toolbox required)
    fs = 44100; % must match soundcard setting

    devices = getAudioDevices(audioDeviceReader);

    device_number = 2; % CHOOSE BEFORE RUN !!!!!!
    device_name = devices{device_number};
    % to check the device number run the following line code in the terminal /
    % command window of matlab and then choose the column number of
    % the USB audio device (sound card):
    % getAudioDevices(audioDeviceReader)
    deviceReader = audioDeviceReader('Device', device_name, 'SampleRate', fs, 'NumChannels', 2, 'SamplesPerFrame', round(Tp*fs));
    disp('Real-time acquisition started');
else
    % pseudo real-time: read wav and feed in chunks
    [data, fs] = audioread(filename);
    % separate channels
    radar_data = data(:,2); % channel2 = radar backscatter data
    sync_data = data(:,1); % channel1 = sync data (square wave)
    disp('Pseudo real-time playback started');
end

%% VIDEO SETUP
saveVideo = true;              % ON(true)/OFF(false) to save/don't save the video
if saveVideo
    v = VideoWriter('FMCW_videos/FMCW_Spectrogram.avi');
    v.FrameRate = 10;           % frames per second
    open(v);
end


%% PLOT INITIALIZATION
N = round(Tp * fs);
fft_size = 4*N;
f_axis = linspace(0, fs/2, fft_size/2);
r_axis = f_axis*c / (2* (stop_freq-start_freq)/Tp);

fig = figure('Color','w', 'Units','pixels', 'Position',[100 100 1074 646], 'Resize','off');
pause(0.1); drawnow;  % let it render once

hImg = imagesc(r_axis, 0, zeros(1,length(r_axis)), [-55 0]);
axis xy;
xlabel('Range (m)');
ylabel('Time (s)');
set(gca,'YDir','reverse');
title('Real-time / Pseudo real-time Range-Time Spectrogram');
colorbar;
xlim([0 50]);   % range axis limit

%% ESC key stop control
stopFlag = false;
set(fig, 'KeyPressFcn', @(src, event) assignin('base','stopFlag', strcmp(event.Key,'escape')));

%% PROCESSING LOOP
frameCount = 0;
t_axis = [];    % store time axis
range_map = [];   % store spectrogram rows
prev_frame = [];

while ~stopFlag
    % ACQUIRE FRAME
    if mode == "real"
        frameStereo = deviceReader();  % Nx2
        if isempty(frameStereo)
            disp("No frame received."); continue;
        end
        frame = frameStereo(:,2);  % take channel 2 (radar signal)
    else
        frameCount = frameCount + 1;
        idx_start = (frameCount-1)*N + 1; %% isnt it same as idx_start = frameCount*N + 1
        idx_end   = frameCount*N;
        if idx_end > length(radar_data), break; end
        frame = radar_data(idx_start:idx_end);
    end

    % PROCESS FRAME
    % --- Mean Subtraction / Clutter Removal ---
    frame = frame - mean(frame); 

    % --- 2-PULSE MTI ---
    if ~isempty(prev_frame)
    mti_frame = frame - prev_frame;
    else
    mti_frame = frame;
    end
    prev_frame = frame;

    % --- FFT PROCESSING ---
    w = hann(N);
    Y = fft(mti_frame .* w, fft_size);
    Y_dB = 20*log10(abs(Y(1:fft_size/2)));
    Y_norm = Y_dB - max(Y_dB);

    % --- UPDATE STORAGE ---
    t_axis(end+1) = frameCount * Tp; 
    range_map(end+1,:) = Y_norm; 

    % --- UPDATE PLOT ---
    set(hImg, 'XData', r_axis, 'YData', t_axis, 'CData', range_map);
    ylim([0 t_axis(end)]); drawnow;

    % --- SAVE VIDEO FRAME ---
    if saveVideo
        try
            frameObj = getframe(fig);
            writeVideo(v, frameObj);
        catch ME
            warning('Skipping frame: %s', ME.message);
            if ~ishandle(fig), break; end
        end
    end


    % --- Pseudo exit condition ---
    if mode == "pseudo" && idx_end >= length(radar_data)
        break;
    end
end


%% === CLEANUP ===
if saveVideo
    close(v);
    disp('Video saved .');
end
