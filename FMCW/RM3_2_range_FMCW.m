%% Range Measurements - FMCW Radar
%% RM3 2. Range Measurements in FMCW mode in FMCW mode

%% Cleaning
clear;
clc;

%% Read data
[data, Fs] = audioread("TEST_FILES18/test1_17_09.wav");

% constants
Tp = 5e-3;            % s
start_freq = 2.405e9; % Hz
stop_freq  = 2.495e9; % Hz
threshold = 0;        % sync threshold
c = 299792458;        % m/s

% separate channels
radar_data = data(:,2); % channel2 = radar backscatter data
sync_data = data(:,1); % channel1 = sync data (square wave)

% %% captured data plot, captured sync plot
% t = (0:length(radar_data)-1)/Fs;
% 
% figure(60);
% subplot(2,1,1);
% plot(t, radar_data);
% title('Radar Backscatter Data');
% xlabel('Time (s)');
% ylabel('Amplitude');
% xlim([2 2.05]);
% grid on;
% 
% subplot(2,1,2);
% plot(t, sync_data);
% title('Sync Signal (Square Wave)');
% xlabel('Time (s)');
% ylabel('Amplitude');
% xlim([2 2.05]);
% grid on;


%% logical sync vector (upchirp when sync > threshold get 1, otherwise 0)
sync_data = sync_data > threshold;

%% samples per chirp
N = floor(Tp*Fs);

%% extract each upchirp using the rising-edge indices
 
% METHOD 2: if the up chirp is not complete, remain data is filled with 0

d_sync = diff([0; sync_data; 0]); %add zeros incase mid up/down chirp for detection
start_idx = find(d_sync == 1);
end_idx = find(d_sync == -1) - 1;

M = numel(start_idx);

%% matrix
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
range_up_CR = matrix_radar_up - mean(matrix_radar_up,1);

% 2-pulse MTI in range-cell columns
% 2-pulse MIT: x[n+1] - x[n] -> produce M-1 rows
radar_up_2MTI = range_up_CR(2:end, :) - range_up_CR(1:end-1,:);

%% IFFT
Y_2MTI = ifft(radar_up_2MTI, 4*N, 2);
Y_2MTI_dB = 20*log10(abs(Y_2MTI));
Y_2MTI_dB_half = Y_2MTI_dB(:, 1:(4*N/2));

%% normalization
max_2MTI = max(Y_2MTI_dB_half(:));
Y_2MTI_dB_half_norm = Y_2MTI_dB_half - max_2MTI; 

%% range axis
B = stop_freq - start_freq;
slope = B/Tp; % Hz/s

% frequency axis for the half-band
df = Fs/(4*N);
f_half = (1:(4*N/2)-1)*df; % 0-based frequency bins (length 4*N/2)
range_axis = (c .* f_half) ./ (2 * slope); % meters
Rmax_check = N * (c/(2*B)) / 2;


%% Time axis (one per matrix)
t_up = chirp_start_times; % M elements
t_2MTI = t_up(1:end-1);  % corresponds to rows of 2MTI (M-1)

%% Range with MS and 2-pulse MTI
figure('Name', 'RM1_FMCW_MS_2MTI');
imagesc(range_axis, t_2MTI, Y_2MTI_dB_half_norm, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
% xlim([0 100]);
xlim([-Inf 30]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] with MS and 2-pulse MTI, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
saveas(gcf, 'FMCW_plots/RM3_2_range_FMCW.png');



%% RM3
%% choose chirps and number of targets
chosen = 10:30;   % change here
%% RM3 with tracking
num_targets = 2;          % track up to 2 targets
range_tracks = nan(M-1, num_targets);

min_sep_bin  = 20;        % minimum bin separation within chirp
max_mov_bin  = 5;        % maximum movement allowed between chirps (bins)

% tracker state
last_bins = nan(1, num_targets);   % last known bin index for each target

for m = 1:M-1
    row_data = Y_2MTI_dB_half_norm(m,:);

    % Step 1: find separated peaks in this chirp
    [pks, locs] = findpeaks(row_data, ...
        'MinPeakDistance', min_sep_bin);

    if isempty(pks)
        continue;
    end

    % sort candidates by strength (descending)
    [~, order] = sort(pks,'descend');
    locs_sorted = locs(order);

    % Step 2: assign to each target ID
    used = false(size(locs_sorted));
    for id = 1:num_targets
        assigned = false;

        if ~isnan(last_bins(id))
            % try to find a candidate close to last known bin
            diffs = abs(locs_sorted - last_bins(id));
            [min_diff, idx] = min(diffs);

            if min_diff <= max_mov_bin
                range_tracks(m,id) = range_axis(locs_sorted(idx));
                last_bins(id) = locs_sorted(idx);
                used(idx) = true;
                assigned = true;
            end
        end

        % if not assigned, pick the strongest unused candidate
        if ~assigned
            idx = find(~used,1,'first');
            if ~isempty(idx)
                range_tracks(m,id) = range_axis(locs_sorted(idx));
                last_bins(id) = locs_sorted(idx);
                used(idx) = true;
            end
        end
    end
end


%% scatterer per chirp
[max_vals, idx_max] = max(Y_2MTI_dB_half_norm, [], 2); 

%% (1) Captured data & sync for chosen chirps
figure('Name','RM3_captured_target');
hold on;
for i = 1:length(chosen)
    m = chosen(i);
    seg_t = (start_idx(m):start_idx(m)+N-1)/Fs;
    seg_radar = radar_data(start_idx(m):start_idx(m)+N-1);
    seg_sync  = sync_data(start_idx(m):start_idx(m)+N-1);

    plot(seg_t, seg_radar,'b');
    plot(seg_t, seg_sync,'r');
end
xlabel('Time (s)'); ylabel('Amplitude');
title(sprintf('Radar Data & Sync (Chirps %d to %d)', chosen(1), chosen(end)));
grid on;
legend('Radar Backscatter','Sync Signal');
saveas(gcf,'FMCW_plots/RM3_captured_chosen.png');

%% (2) Overlay radar + sync for chosen chirps
figure('Name','RM3_overlay');
for i = 1:length(chosen)
    m = chosen(i);
    seg_t = (start_idx(m):start_idx(m)+N-1)/Fs;
    seg_radar = radar_data(start_idx(m):start_idx(m)+N-1);
    seg_sync  = sync_data(start_idx(m):start_idx(m)+N-1);

    yyaxis left
    plot(seg_t, seg_radar, 'b'); ylabel('Radar Data');
    yyaxis right
    plot(seg_t, seg_sync, 'r'); ylim([0 1]); ylabel('Sync Signal');
end
xlabel('Time (s)');
title(sprintf('Overlay Radar & Sync (Chirps %d to %d)', chosen(1), chosen(end)));
grid on;
saveas(gcf,'FMCW_plots/RM3_overlay_chosen.png');

%% (3) Range vs Time (all strongest scatterers track)
figure('Name','RM3_range_vs_time');
plot(t_2MTI, range_tracks(:,1),'ro'); hold on;
if num_targets >= 2
    plot(t_2MTI, range_tracks(:,2),'cx');
end
% plot(t_2MTI(chosen), range_track(chosen),'ro','MarkerSize',4,'LineWidth',1.5);
xlabel('Time (s)'); ylabel('Range (m)');
title(sprintf('Range vs Time (%d Strongest Scatterers)', num_targets));
legend_strings = arrayfun(@(k) sprintf('Target %d',k),1:num_targets,'UniformOutput',false);
legend(legend_strings);
grid on;
saveas(gcf,'FMCW_plots/RM3_range_vs_time_multi.png');

%% (4) Range-Time Spectrogram + strongest target
figure('Name','RM3_spectrogram');
imagesc(t_2MTI, range_axis, Y_2MTI_dB_half_norm.',[-50 0]); 
axis xy; colorbar;
xlabel('Time (s)'); ylabel('Range (m)');
ylim([-Inf 30]);
title(sprintf('Range-Time Spectrogram with %d Targets', num_targets));
hold on;
plot(t_2MTI, range_tracks(:,1),'ro');

% plot(t_2MTI(chosen), range_track(chosen),'ro','MarkerSize',4,'LineWidth',1.5);
if num_targets >= 2
    plot(t_2MTI, range_tracks(:,2),'cx');
end
legend(legend_strings);
saveas(gcf,'FMCW_plots/RM3_range_vs_time_multi.png');
