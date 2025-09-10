%% Range Measurements - FMCW Radar
% 
% Range Algorithm development for FMCW Radar.
% Variation in results with and without MS clutter rejection.
% Variation in results utilizing more, less or no zero padding.
% Results for both 2-pulse and 3-pulse MTI clutter rejection.
% 
% Group 7
% 10/09/2025


%% Cleaning
clear;
clc;

%% Read data
[data, Fs] = audioread("Range_Test_File2.m4a");  % user file

% Constants
Tp = 5e-3;            % s
start_freq = 2.405e9; % Hz
stop_freq  = 2.495e9; % Hz
threshold = 0;        % sync threshold
c = 299792458;        % m/s

% Separate channels
radar_data = data(:,1); % channel 1 = radar backscatter data
sync_data  = data(:,2); % channel 2 = sync data (square wave)

%% Check the data
%t = (0:length(radar_data)-1)/Fs;

%figure;
%subplot(2,1,1);
%plot(t, radar_data);
%title('Radar Backscatter Data');
%xlabel('Time (s)');
%ylabel('Amplitude');
%xlim([2 2.05]);
%grid on;

%subplot(2,1,2);
%plot(t, sync_data);
%title('Sync Signal (Square Wave)');
%xlabel('Time (s)');
%ylabel('Amplitude');
%xlim([2 2.05]);
%grid on;

%% Logical sync vector (up-chirp when sync > threshold get 1, otherwise 0)

sync_up = sync_data > threshold;

% Imperfect sync wave: 

% t = (0:length(radar_data)-1)/Fs;
% start_idx = 3265043;
% end_idx   = 3265090;
% figure(30);
% plot(t(start_idx:end_idx), sync_data(start_idx:end_idx));
% title('sync wave between selected indexes: 3265043 and 3265090');
% xlabel('Time (s)');
% ylabel('sync wave');
% grid on;

%% Samples per chirp
N = Tp * Fs;

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

% Optional: to process the downchirp it is the same process, with falling edges
% falling = find(diff([sync_up; 0]) == -1);

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
mean_range_cell = mean(matrix_radar_up, 1);
range_up_CR = matrix_radar_up - mean_range_cell;

%% 2-pulse and 3-pulse MTI (in range-cell columns)
% 2-pulse MTI: x[n+1] - x[n] -> produce M-1 rows
radar_up_2MTI = range_up_CR(2:end, :) - range_up_CR(1:end-1, :);

% 3-pulse: x[n+2] - 2*x[n+1] + x[n] -> produces M-2 rows
radar_up_3MTI = range_up_CR(3:end, :) - 2*range_up_CR(2:end-1, :) + range_up_CR(1:end-2, :);

%% IFFT
Y_raw = ifft(matrix_radar_up, 4*N, 2);
Y_CR  = ifft(range_up_CR, 4*N, 2);
Y_2MTI = ifft(radar_up_2MTI, 4*N, 2);
Y_3MTI = ifft(radar_up_3MTI, 4*N, 2);

% Less zero padding and no zero padding / more zero padding
Y_raw_lzp = ifft(matrix_radar_up, 2*N, 2);
Y_CR_lzp = ifft(range_up_CR, 2*N, 2);
Y_2MTI_lzp = ifft(radar_up_2MTI, 2*N, 2);
Y_3MTI_lzp = ifft(radar_up_3MTI, 2*N, 2);

Y_raw_nzp = ifft(matrix_radar_up, [], 2);
Y_CR_nzp = ifft(range_up_CR, [], 2);
Y_2MTI_nzp = ifft(radar_up_2MTI, [], 2);
Y_3MTI_nzp = ifft(radar_up_3MTI, [], 2);

Y_CR_mzp  = ifft(range_up_CR, 8*N, 2);
% end

Y_raw_dB = 20*log10(abs(Y_raw));
Y_CR_dB  = 20*log10(abs(Y_CR));
Y_2MTI_dB = 20*log10(abs(Y_2MTI));
Y_3MTI_dB = 20*log10(abs(Y_3MTI));

% Less zero padding and no zero padding / more zero padding
Y_raw_dB_lzp = 20*log10(abs(Y_raw_lzp));
Y_CR_dB_lzp = 20*log10(abs(Y_CR_lzp));
Y_2MTI_dB_lzp = 20*log10(abs(Y_2MTI_lzp));
Y_3MTI_dB_lzp = 20*log10(abs(Y_3MTI_lzp));

Y_raw_dB_nzp = 20*log10(abs(Y_raw_nzp));
Y_CR_dB_nzp = 20*log10(abs(Y_CR_nzp));
Y_2MTI_dB_nzp = 20*log10(abs(Y_2MTI_nzp));
Y_3MTI_dB_nzp = 20*log10(abs(Y_3MTI_nzp));

Y_CR_dB_mzp = 20*log10(abs(Y_CR_mzp));
% end


Y_raw_dB_half = Y_raw_dB(:, 1:(4*N/2));
Y_CR_dB_half  = Y_CR_dB(:, 1:(4*N/2));
Y_2MTI_dB_half = Y_2MTI_dB(:, 1:(4*N/2));
Y_3MTI_dB_half = Y_3MTI_dB(:, 1:(4*N/2));

% Less zero padding and no zero padding / more zero padding
Y_raw_dB_half_lzp = Y_raw_dB_lzp(:, 1:(2*N/2));
Y_CR_dB_half_lzp  = Y_CR_dB_lzp(:, 1:(2*N/2));
Y_2MTI_dB_half_lzp = Y_2MTI_dB_lzp(:, 1:(2*N/2));
Y_3MTI_dB_half_lzp = Y_3MTI_dB_lzp(:, 1:(2*N/2));

Y_raw_dB_half_nzp = Y_raw_dB_nzp(:, 1:(N/2));
Y_CR_dB_half_nzp  = Y_CR_dB_nzp(:, 1:(N/2));
Y_2MTI_dB_half_nzp = Y_2MTI_dB_nzp(:, 1:(N/2));
Y_3MTI_dB_half_nzp = Y_3MTI_dB_nzp(:, 1:(N/2));

Y_CR_dB_half_mzp  = Y_CR_dB_mzp(:, 1:(8*N/2));
% end

%% Normalization
max_raw  = max(Y_raw_dB_half(:));
max_CR   = max(Y_CR_dB_half(:));
max_2MTI = max(Y_2MTI_dB_half(:));
max_3MTI = max(Y_3MTI_dB_half(:));

% Less zero padding and no zero padding / more zero padding
max_raw_lzp  = max(Y_raw_dB_half_lzp(:));
max_CR_lzp   = max(Y_CR_dB_half_lzp(:));
max_2MTI_lzp = max(Y_2MTI_dB_half_lzp(:));
max_3MTI_lzp = max(Y_3MTI_dB_half_lzp(:));

max_raw_nzp  = max(Y_raw_dB_half_nzp(:));
max_CR_nzp   = max(Y_CR_dB_half_nzp(:));
max_2MTI_nzp = max(Y_2MTI_dB_half_nzp(:));
max_3MTI_nzp = max(Y_3MTI_dB_half_nzp(:));

max_CR_mzp   = max(Y_CR_dB_half_mzp(:));
% end


Yup_raw_dB_norm  = Y_raw_dB_half  - max_raw;
Yup_CR_dB_norm   = Y_CR_dB_half   - max_CR;
Yup_2MTI_dB_norm = Y_2MTI_dB_half - max_2MTI;
Yup_3MTI_dB_norm = Y_3MTI_dB_half - max_3MTI;

% Less zero padding and no zero padding / more zero padding
Yup_raw_dB_norm_lzp  = Y_raw_dB_half_lzp  - max_raw_lzp;
Yup_CR_dB_norm_lzp   = Y_CR_dB_half_lzp   - max_CR_lzp;
Yup_2MTI_dB_norm_lzp = Y_2MTI_dB_half_lzp - max_2MTI_lzp;
Yup_3MTI_dB_norm_lzp = Y_3MTI_dB_half_lzp - max_3MTI_lzp;

Yup_raw_dB_norm_nzp  = Y_raw_dB_half_nzp  - max_raw_nzp;
Yup_CR_dB_norm_nzp   = Y_CR_dB_half_nzp   - max_CR_nzp;
Yup_2MTI_dB_norm_nzp = Y_2MTI_dB_half_nzp - max_2MTI_nzp;
Yup_3MTI_dB_norm_nzp = Y_3MTI_dB_half_nzp - max_3MTI_nzp;

Yup_CR_dB_norm_mzp   = Y_CR_dB_half_mzp   - max_CR_mzp;
% end

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

%% Plotting (imagesc expects x-axis = range, y-axis = time)
figure(1);
imagesc(range_axis, t_up, Yup_raw_dB_norm, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] without MS and no MTI, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
% saveas(gcf, 'Raw_up.png');

figure(2);
imagesc(range_axis, t_up, Yup_CR_dB_norm, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] with MS and no MTI, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
% saveas(gcf, 'MS_up.png');

figure(3);
imagesc(range_axis, t_2MTI, Yup_2MTI_dB_norm, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] with MS and 2-pulse MTI, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
% saveas(gcf, 'MS_2MTI_up.png');

figure(4);
imagesc(range_axis, t_3MTI, Yup_3MTI_dB_norm, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] with MS and 3-pulse MTI, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
% saveas(gcf, 'MS_3MTI_up.png');


% Less zero padding
df = Fs / (2*N);
f_half = ((1:(2*N/2)) - 1) * df;   % 0-based frequency bins (length 2*N/2)
range_axis = (c .* f_half) ./ (2 * slope);   % meters

figure(5);
imagesc(range_axis, t_up, Yup_raw_dB_norm_lzp, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] without MS and no MTI, less zero padding, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
% saveas(gcf, 'Raw_up_lzp.png');

figure(6);
imagesc(range_axis, t_up, Yup_CR_dB_norm_lzp, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] with MS and no MTI, less zero padding, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
% saveas(gcf, 'MS_up_lzp.png');

figure(7);
imagesc(range_axis, t_2MTI, Yup_2MTI_dB_norm_lzp, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] with MS and 2-pulse MTI, less zero padding, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
% saveas(gcf, 'MS_2MTI_up_lzp.png');

figure(8);
imagesc(range_axis, t_3MTI, Yup_3MTI_dB_norm_lzp, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] with MS and 3-pulse MTI, less zero padding, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
% saveas(gcf, 'MS_3MTI_up_lzp.png');


% No zero padding
df = Fs / (N);
f_half = ((1:(N/2)) - 1) * df;   % 0-based frequency bins (length N/2)
range_axis = (c .* f_half) ./ (2 * slope);   % meters

figure(9);
imagesc(range_axis, t_up, Yup_raw_dB_norm_nzp, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] without MS and no MTI, no zero padding, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
% saveas(gcf, 'Raw_up_nzp.png');

figure(10);
imagesc(range_axis, t_up, Yup_CR_dB_norm_nzp, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] with MS and no MTI, no zero padding, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
% saveas(gcf, 'MS_up_nzp.png');

figure(11);
imagesc(range_axis, t_2MTI, Yup_2MTI_dB_norm_nzp, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] with MS and 2-pulse MTI, no zero padding, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
% saveas(gcf, 'MS_2MTI_up_nzp.png');

figure(12);
imagesc(range_axis, t_3MTI, Yup_3MTI_dB_norm_nzp, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] with MS and 3-pulse MTI, no zero padding, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
% saveas(gcf, 'MS_3MTI_up_nzp.png');


% More zero padding
df = Fs / (8*N);
f_half = ((1:(8*N/2)) - 1) * df;   % 0-based frequency bins (length N/2)
range_axis = (c .* f_half) ./ (2 * slope);   % meters

figure(13);
imagesc(range_axis, t_3MTI, Yup_CR_dB_norm_mzp, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] with MS and no MTI, more zero padding, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
% saveas(gcf, 'MS_3MTI_up_mzp.png');

