%% Range Measurements - FMCW Radar

%% Cleaning
clear;
clc;

%% Read data
[data, Fs] = audioread("test1_17_09.wav");  % user file

% Constants
Tp = 5e-3;            % s
start_freq = 2.405e9; % Hz
stop_freq  = 2.495e9; % Hz
threshold = 0;        % sync threshold
c = 299792458;        % m/s

% Separate channels
radar_data = data(:,2); % channel 2 = radar backscatter data
sync_data  = data(:,1); % channel 1 = sync data (square wave)

%% Check the data
t = (0:length(radar_data)-1)/Fs;

figure(60);
subplot(2,1,1);
plot(t, radar_data);
title('Radar Backscatter Data');
xlabel('Time (s)');
ylabel('Amplitude');
% xlim([2 2.05]);
grid on;

subplot(2,1,2);
plot(t, sync_data);
title('Sync Signal (Square Wave)');
xlabel('Time (s)');
ylabel('Amplitude');
% xlim([2 2.05]);
grid on;

%% Logical sync vector (up-chirp when sync > threshold get 1, otherwise 0)

sync_up = sync_data > threshold;


%% Samples per chirp
N = floor(Tp * Fs);

%% Extract each up-chirp using the rising-edge indices

% METHOD 1: get N points even if up chirp is not complete

% Find rising edges: indices where sync goes 0 -> 1
% rising = find(diff([0; sync_up]) == 1); % indexes of the first sample where up starts
% M = length(rising);
% 
% matrix_radar_up = zeros(M, N);
% chirp_start_times = zeros(M,1);
% 
% for m = 1:M
%     idx = rising(m);
%     matrix_radar_up(m, :) = radar_data(idx : idx + N - 1).'; % row = one chirp
%     chirp_start_times(m) = (idx - 1) / Fs; % seconds
% end
% 
% % Optional: to process the downchirp it is the same process, with falling edges
% % falling = find(diff([sync_up; 0]) == -1);
% 
% METHOD 2: if the up chirp is not complete, remain data is filled with 0

d_sync = diff([0; sync_up; 0]); %add zeros incase mid up/down chirp for detection
start_idx = find(d_sync == 1);
end_idx = find(d_sync == -1) - 1;

M = numel(start_idx);

matrix_radar_up = zeros(M, N);
chirp_start_times = zeros(M,1);

for m = 1:M
    if end_idx(m) - start_idx(m) + 1 >= N % complete up-chirp
        matrix_radar_up(m,:) = radar_data(start_idx(m) : start_idx(m) + N - 1);
    else
        % if a chirp is shorter than expected -> 0-pad
        seg = radar_data(start_idx(m):end_idx(m));
        matrix_radar_up(m,1:numel(seg)) = seg;
    end
    chirp_start_times(m) = (start_idx(m) - 1) / Fs; % seconds
end

%% MS clutter rejection: subtract mean across chirps for each range cell (column)
range_up_CR = matrix_radar_up - mean(matrix_radar_up, 1);

% 2-pulse and 3-pulse MTI (in range-cell columns)
% 2-pulse MTI: x[n+1] - x[n] -> produce M-1 rows
radar_up_2MTI = range_up_CR(2:end, :) - range_up_CR(1:end-1, :);

% 3-pulse: x[n+2] - 2*x[n+1] + x[n] -> produces M-2 rows
radar_up_3MTI = range_up_CR(3:end, :) - 2*range_up_CR(2:end-1, :) + range_up_CR(1:end-2, :);

%% IFFT
Y_raw = ifft(matrix_radar_up, 4*N, 2);
Y_CR  = ifft(range_up_CR, 4*N, 2);
Y_2MTI = ifft(radar_up_2MTI, 4*N, 2);
Y_2MTI_fft = fft(radar_up_2MTI, 4*N, 2);

Y_3MTI = ifft(radar_up_3MTI, 4*N, 2);

Y_raw_dB = 20*log10(abs(Y_raw));
Y_CR_dB  = 20*log10(abs(Y_CR));
Y_2MTI_dB = 20*log10(abs(Y_2MTI));
Y_2MTI_dB_fft = 20*log10(abs(Y_2MTI_fft));

Y_3MTI_dB = 20*log10(abs(Y_3MTI));

Y_raw_dB_half = Y_raw_dB(:, 1:(4*N/2));
Y_CR_dB_half  = Y_CR_dB(:, 1:(4*N/2));
Y_2MTI_dB_half = Y_2MTI_dB(:, 1:(4*N/2));
Y_2MTI_dB_half_fft = Y_2MTI_dB_fft(:, 1:(4*N/2));

Y_3MTI_dB_half = Y_3MTI_dB(:, 1:(4*N/2));


%% Normalization
max_raw  = max(Y_raw_dB_half(:));
max_CR   = max(Y_CR_dB_half(:));
max_2MTI = max(Y_2MTI_dB_half(:));
max_2MTI_fft = max(Y_2MTI_dB_half_fft(:));
max_3MTI = max(Y_3MTI_dB_half(:));

Yup_raw_dB_norm  = Y_raw_dB_half  - max_raw;
Yup_CR_dB_norm   = Y_CR_dB_half   - max_CR;
Yup_2MTI_dB_norm = Y_2MTI_dB_half - max_2MTI;
Yup_2MTI_dB_norm_fft = Y_2MTI_dB_half_fft - max_2MTI_fft;
Yup_3MTI_dB_norm = Y_3MTI_dB_half - max_3MTI;


%% Range axis
B = stop_freq - start_freq;
slope = B / Tp;                  % Hz/s

% frequency axis for the half-band
df = Fs / (4*N);
f_half = (1:(4*N/2) - 1) * df; % 0-based frequency bins (length 4*N/2)
range_axis = (c .* f_half) ./ (2 * slope); % meters

Rmax_check = N * (c/(2*B)) / 2;

%% Time axes (one per matrix)
t_up = chirp_start_times; % M elements
t_2MTI = t_up(1:end-1); % corresponds to rows of 2MTI (M-1)
t_3MTI = t_up(1:end-2); % corresponds to rows of 3MTI (M-2)

% %% Plotting (imagesc expects x-axis = range, y-axis = time)
% figure(1);
% imagesc(range_axis, t_up, Yup_raw_dB_norm, [-50 0]);
% axis xy;
% colorbar;
% set(gca,'YDir','reverse');
% % xlim([0 100]);
% xlim([-Inf 30]);
% xlabel('Range (m)');
% ylabel('Time (s)');
% title('Range [upc] without MS and no MTI, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
% saveas(gcf, 'FMCW_plots/Raw_up.png');

% figure(2);
% imagesc(range_axis, t_up, Yup_CR_dB_norm, [-50 0]);
% axis xy;
% colorbar;
% set(gca,'YDir','reverse');
% % xlim([0 100]);
% xlim([-Inf 30]);
% xlabel('Range (m)');
% ylabel('Time (s)');
% title('Range [upc] with MS and no MTI, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
% saveas(gcf, 'FMCW_plots/MS_up.png');
% 
figure('Name', '3_ifft');
imagesc(range_axis, t_2MTI, Yup_2MTI_dB_norm, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
% xlim([0 100]);
xlim([-Inf 30]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] with MS and 2-pulse MTI, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
saveas(gcf, 'FMCW_plots/MS_2MTI_up_ifft.png');

figure('Name', '3_fft');
imagesc(range_axis, t_2MTI, Yup_2MTI_dB_norm_fft, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
% xlim([0 100]);
xlim([-Inf 30]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] with MS and 2-pulse MTI, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
saveas(gcf, 'FMCW_plots/MS_2MTI_up_fft.png');


% figure(4);
% imagesc(range_axis, t_3MTI, Yup_3MTI_dB_norm, [-50 0]);
% axis xy;
% colorbar;
% set(gca,'YDir','reverse');
% % xlim([0 100]);
% xlim([-Inf 30]);
% xlabel('Range (m)');
% ylabel('Time (s)');
% title('Range [upc] with MS and 3-pulse MTI, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
% saveas(gcf, 'FMCW_plots/MS_3MTI_up.png');
% 
% 
