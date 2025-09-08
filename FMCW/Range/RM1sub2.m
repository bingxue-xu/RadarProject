%% Cleaning
clear;
clc;

%% Get datas
[data, Fs] = audioread("Range_Test_File2.wav");

Tp = 5e-3; % ms
start_freq = 2.405e9; % Hz
stop_freq = 2.495e9; % Hz
threshold = 0; % sync threshold
c = 299792458; % m/s

% Separate channels
radar_data = data(:,1);   % channel 1 = radar backscatter data
sync_data  = data(:,2);   % assume channel 2 = sync data (square wave)

%% Just to check the data
% Create time axis
%t = (0:length(radar_data)-1)/Fs;

% Plot both signals
%figure;

%subplot(2,1,1);
%plot(t, radar_data);
%title('Radar Backscatter Data');
%xlabel('Time (s)');
%ylabel('Amplitude');
%grid on;

%subplot(2,1,2);
%plot(t, sync_data);
%title('Sync Signal (Square Wave)');
%xlabel('Time (s)');
%ylabel('Amplitude');
%grid on;


%% Extract chirps

sync_pulse_up = sync_data > threshold;   % up-chirp sync
%sync_pulse_down  = ~sync_pulse_up;       % down-chirp sync

radar_up_data = radar_data(sync_pulse_up == 1);
%radar_down_data = radar_data(sync_pulse_down == 1);

N_upchirps = length(radar_up_data);
%N_downchirps = length(radar_down_data);

%% Processing only up-chirp

N = Tp*Fs;
M = ceil(length(radar_up_data)/N);

matrix_radar_up = zeros(M,N);

for i = 0:M-1
    for j = 1:N
        if i*N+j <= length(radar_up_data)
            matrix_radar_up(i+1,j) = radar_up_data(i*N+j,1); % Inserting element in matrix
        end
    end
end


%% Perform MS clutter rejection
mean_range_cell = mean(matrix_radar_up, 1); % mean of each column

range_up_CR = matrix_radar_up - mean_range_cell;

radar_up_2MTI = range_up_CR(2:size(range_up_CR,1),:) - range_up_CR(1:size(range_up_CR,1)-1,:);

radar_up_3MTI = range_up_CR(3:size(range_up_CR,1),:) - range_up_CR(2:size(range_up_CR,1)-1,:) - ( range_up_CR(2:size(range_up_CR,1)-1,:) - range_up_CR(1:size(range_up_CR,1)-2,:) );

%% IFFT
Y_raw = ifft(matrix_radar_up, 4*N, 2);
Yup_CR = ifft(range_up_CR, 4*N, 2);
Yup_2MTI = ifft(radar_up_2MTI, 4*N, 2);
Yup_3MTI = ifft(radar_up_3MTI, 4*N, 2);

Yup_raw_dB = 20*log10(abs(Y_raw));
Yup_CR_dB = 20*log10(abs(Yup_CR));
Yup_2MTI_dB = 20*log10(abs(Yup_2MTI));
Yup_3MTI_dB = 20*log10(abs(Yup_3MTI));

% Calculating Fmax value
Fmax = Fs/2;

Yup_raw_dB_half = Yup_raw_dB(:,1:2*N);
Yup_CR_dB_half = Yup_CR_dB(:,1:2*N);
Yup_2MTI_dB_half = Yup_2MTI_dB(:,1:2*N);
Yup_3MTI_dB_half = Yup_3MTI_dB(:,1:2*N);


%% Normalization 1
% perform normalization on this half matrix by subtracting 
% the maximum of the matrix from each value of the matrix.

% Find the maximum for each row, and then the maximum among these values
max_raw = max(Yup_raw_dB_half(:));
max_CR = max(Yup_CR_dB_half(:));
max_2MTI = max(Yup_2MTI_dB_half(:));
max_3MTI = max(Yup_3MTI_dB_half(:));

% Subtracting the maximum value
Yup_raw_dB_norm = Yup_raw_dB_half - max_raw;
Yup_CR_dB_norm = Yup_CR_dB_half - max_CR;
Yup_2MTI_dB_norm = Yup_2MTI_dB_half - max_2MTI;
Yup_3MTI_dB_norm = Yup_3MTI_dB_half - max_3MTI;


%% Axis
delta_f = stop_freq - start_freq;
delta_R = c/(2*delta_f);
Rmax = N*delta_R/2;

f_axis = linspace(0, Fmax, 2*N);
x_axis = linspace(0, Rmax, 2*N);
t_axis = linspace(1, Tp*M, M);

figure(1)
imagesc(x_axis, t_axis, Yup_raw_dB_norm, [-50 0])
colorbar;
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range without MS and no MTI, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');

figure(2)
imagesc(x_axis, t_axis, Yup_CR_dB_norm, [-50 0])
colorbar;
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range with MS and no MTI, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');

figure(3)
imagesc(x_axis, t_axis, Yup_2MTI_dB_norm, [-50 0])
colorbar;
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range with MS and 2-pulse MTI, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');

figure(4)
imagesc(x_axis, t_axis, Yup_3MTI_dB_norm, [-50 0])
colorbar;
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range with MS and 3-pulse MTI, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');





