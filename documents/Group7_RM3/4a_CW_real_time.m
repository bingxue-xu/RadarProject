%% Real-time Measurements - CW Mode
% 
% Group 7
% Review Meeting 3
% 23/09/2025
% 

%% READ:
% 
% In pseudo real-time mode with file, wait for execution to finish so as 
% not to corrupt the video
% 
% In real-time mode:
% 1. Identify the correct input device to acquire the data
% 2. To finish acquisition/video, simply press the ‘ESC’ key to exit the 
% program correctly and avoid corrupting the video
% 

%% CW Doppler Radar – Real-time / Pseudo real-time velocity measurements
clear; clc; close all;

%% USER SETTINGS
mode = "real";   % choose: "real" or "pseudo"
filename = "Test_files/CW_2targets_t2.wav"; % only used in pseudo mode
Tp = 0.1;                 % processing window (s)
fc = 2.445e9;              % radar center frequency (Hz)
c = 3e8;                   % speed of light

%% ACQUISITION SETUP
if mode == "real"
    % live from soundcard (Audio Toolbox required)
    fs = 44100; % must match soundcard setting

    devices = getAudioDevices(audioDeviceReader);

    device_number = 3; % CHOOSE BEFORE RUN!
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
    data = data(:,2);  % assuming radar channel is right channel
    N = round(Tp*fs);
    M = floor(length(data)/N);
    disp('Pseudo real-time playback started');
end

%% VIDEO SETUP
saveVideo = true;              % ON(true)/OFF(false) to save/don't save the video
if saveVideo
    v = VideoWriter('CW_videos/CW_Spectrogram.avi');
    v.FrameRate = 10;           % frames per second
    open(v);
end

%% FFT PARAMETERS
N = round(Tp*fs);
fft_size = 4*N;
f_axis = linspace(0, fs/2, fft_size/2);
v_axis = f_axis*c/(2*fc);

%% PLOT INITIALIZATION
fig = figure('Color','w');
% fig.WindowState = 'maximized';
hImg = imagesc(v_axis, 0, zeros(1,length(v_axis)), [-55 0]); % empty spectrogram
axis xy;
xlabel('Velocity (m/s)');
ylabel('Time (s)');
set(gca,'YDir','reverse');
title('Real-time / Pseudo real-time Velocity-Time Spectrogram');
colorbar;
xlim([0 10]);   % velocity axis limit

%% ESC key stop control
stopFlag = false;
set(fig, 'KeyPressFcn', @(src, event) assignin('base','stopFlag', strcmp(event.Key,'escape')));


%% PROCESSING LOOP
frameCount = 0;
t_axis = [];    % store time axis
vel_map = [];   % store spectrogram rows

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
        idx_start = (frameCount-1)*N + 1;
        idx_end   = frameCount*N;
        if idx_end > length(data), break; end
        frame = data(idx_start:idx_end);
    end

    % PROCESS FRAME
    frame = frame - mean(frame); % clutter removal
    % frame = frame - mean(frame,1);

    Y = fft(frame, fft_size);
    Y_dB = 20*log10(abs(Y(1:fft_size/2)));

    % Normalize (optional: choose one)
    Y_norm = Y_dB - max(Y_dB);

    % UPDATE DATA STRUCTURES
    frameCount = frameCount + 1;
    % This flags below "%#ok" are used just to eliminate the warnings - if you remove the flags, 
    % warnings appear in the lines (vector that keeps growing at every iteration)
    t_axis(end+1) = frameCount*Tp; %#ok
    vel_map(end+1,:) = Y_norm; %#ok

    % UPDATE PLOT
    set(hImg, 'XData', v_axis, 'YData', t_axis, 'CData', vel_map);
    ylim([-Inf t_axis(end)]);
    
    drawnow;

    % SAVE FRAME TO VIDEO
    if saveVideo
        frameObj = getframe(fig);
        writeVideo(v, frameObj);
    end

    % BREAK CONDITION (pseudo only)
    if mode == "pseudo" && frameCount >= M
        break;
    end
end


if saveVideo
    close(v);
    disp('Video saved successfully.');
end


