%% Range Measurements - FMCW Radar

%% Cleaning
clear;
clc;
close all;

%% Read data
[data, Fs] = audioread("./TEST_FILES/FMCW_19_t1_3.wav");  % user file

% Constants
Tp = 5e-3;            % s
start_freq = 2.4e9; % Hz
stop_freq  = 2.5e9; % Hz
threshold = 0;        % sync threshold
c = 299792458;        % m/s

% Separate channels
radar_data = data(:,2); % channel 1 = radar backscatter data
sync_data  = data(:,1); % channel 2 = sync data (square wave)
shift_amount = floor(Tp*Fs*0.5);
sync_data = [sync_data(shift_amount+1:end); zeros(shift_amount,1)];
%% Check the data
t = (0:length(radar_data)-1)/Fs;

figure(1);
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
range_CR = matrix_radar_up - mean(matrix_radar_up, 1);

%% 3-pulse MTI (in range-cell columns)
% 3-pulse: x[n+2] - 2*x[n+1] + x[n] -> produces M-2 rows
radar_3MTI = range_CR(3:end, :) - 2*range_CR(2:end-1, :) + range_CR(1:end-2, :);

%% IFFT
Y_3MTI = abs(ifft(radar_3MTI, 4*N, 2));

Y_3MTI_dB = 20*log10(Y_3MTI);

Y_3MTI_dB_half = Y_3MTI_dB(:, 1:2*N);

%% Normalization
max_3MTI = max(Y_3MTI_dB_half(:));

Y_3MTI_dB_norm = Y_3MTI_dB_half - max_3MTI;

%% Range axis
B = stop_freq - start_freq;
slope = B / Tp;                  % Hz/s

% frequency axis for the half-band
df = Fs / (4*N);
f_bins = (0:2*N - 1) * df; % 0-based frequency bins (length 4*N/2)
Rmax = N*c/(4*B);
range_axis = linspace(0,Rmax,2*N); %meters

%% Time axes (one per matrix)
t_up = chirp_start_times; % M elements
t_3MTI = t_up(1:end-2); % corresponds to rows of 3MTI (M-2)

%% Plotting (imagesc expects x-axis = range, y-axis = time)

figure(2);
imagesc(range_axis, t_3MTI, Y_3MTI_dB_norm, [-50 0]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
xlim([-Inf 40]);
xlabel('Range (m)');
ylabel('Time (s)');
title('Range [upc] with MS and 3-pulse MTI, Tp=5ms, fstart=2.405GHz, fstop=2.495GHz');
saveas(gcf, 'FMCW_plots/MS_3MTI_up.png');


%% --- Compute beat frequency (strongest peak) ---

% Ignore DC bin (set to zero) to avoid strong DC bias
Y_3MTI_dB_norm(:,1) = -Inf;
% find index of maximum magnitude for each chirp
[vals, idx_max] = max(Y_3MTI_dB_norm, [], 2);   % M x 1 index vector
idx_max(1) = 1; % Throw the first one

%% Filter results by removing weak signals
for i = 2:length(idx_max)
    if vals(i) < -25 %|| abs(idx_max(i)-idx_max(i-1)) > 30
        idx_max(i) = idx_max(i-1);
    end
end

%% Evaluate range with just up-chirp or just down-chirp
% convert to beat frequency (Hz)
f_beat = f_bins(idx_max).';            % 1 x M (transpose to row if you prefer)

range_from_beat = (c .* f_beat) ./ (4*B/(2*Tp));  % m

% Times for chirps
t_chirp = chirp_start_times;  % M x 1 as in your code
%% Plot range with just up-chirp or just down-chirp
% Plot beat frequency vs time and range vs time
figure(4);
subplot(1,2,1);
plot(t_chirp(1:M-2), f_beat, '.-');
xlabel('Time (s)'); ylabel('Beat frequency (Hz)');
title('Beat frequency per up-chirp (strongest peak)');
grid on;

subplot(1,2,2);
plot(t_chirp(1:M-2), range_from_beat, '.-');
xlabel('Time (s)'); ylabel('Range (m)');
title('Estimated range from up beat frequency');
grid on;


saveas(gcf, 'FMCW_plots/Range_and_bfreq.png');


%% --- Velocity estimation ---
period = 0.3; %s
points = period/(2*Tp);
samples = ceil(length(range_from_beat)/points);

vel_est = zeros(1,samples);

for i = 1:samples
    if i*points > length(range_from_beat)
        vel_est(i) = (range_from_beat(end) - range_from_beat((i-1)*points + 1))/period;
    else
        vel_est(i) = (range_from_beat(i*points) - range_from_beat((i-1)*points + 1))/period;
    end
end

%% Filtering by removing too fast variation between consecutive elements
% for i = 2:length(range_est)
%     if abs(range_est(i)-range_est(i-1)) > 2
%         range_est(i) = range_est(i-1);
%     end
% end
% 
% for i = 2:length(vel_est)
%     if abs(vel_est(i)-vel_est(i-1)) > 3
%         vel_est(i) = vel_est(i-1);
%     end
% end

%% Plot range and velocity using beat and doppler frequencies

% Interpolation with fifth order polynomial
t_vel = linspace(0,max(t_3MTI),length(vel_est));
order = 4;
pr = polyfit(t_3MTI, range_from_beat, order);
pv = polyfit(t_vel, vel_est, order);
time_p = linspace(min(t_3MTI), max(t_3MTI), 5000);
range_smooth = polyval(pr, time_p);
vel_smooth = polyval(pv, time_p);

figure(6);
plot(t_3MTI, range_from_beat, '.-');
hold on
plot(time_p, range_smooth, '-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Range (m)');
title('Estimated range from beat frequency');
grid on;
saveas(gcf, 'FMCW_plots/range.png');

figure(7);
plot(t_vel, vel_est, '.-');
hold on
plot(time_p, vel_smooth, '-', 'LineWidth', 2);
hold off
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Estimated velocity from beat frequency');
grid on;
saveas(gcf, 'FMCW_plots/velocity.png');
