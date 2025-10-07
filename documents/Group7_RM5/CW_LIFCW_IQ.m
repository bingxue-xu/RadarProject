%% Velocity Measurements - CW Mode - Single and Multiple Targets
% 
% Group 7
% Review Meeting 5
% 07/10/2025
% 

%% Cleaning
clear;
clc;
close all;

%% Read data
% Read both I and Q data
[data_I, Fs] = audioread("TESTS_02_10/LIF_CW_1target/low_if_cw_re.wav");
[data_Q, ~] = audioread("TESTS_02_10/LIF_CW_1target/low_if_cw_im.wav");

% Combine into complex baseband
data = data_I + 1j*data_Q;

% data = data(0.1*Fs:end);

Tp = 0.1;
fc = 5.8e9;
c = 3e8;
N = Tp*Fs;

M = floor(length(data)/N);

%% Plot captured data
t = (0:length(data)-1)/Fs;  % time axis in seconds

figure(1);
plot(t, data);
xlabel('Time (s)');
ylabel('Amplitude');
title('Captured Radar Signal');
grid on;
% xlim([5 5.05]);
saveas(gcf, 'LOW_IF_CW_1target_plots_IQ/Captured_Signal_total.png');

figure(2);
plot(t, data);
xlabel('Time (s)');
ylabel('Amplitude');
title('Captured Radar Signal');
xlim([5 5.05]);
grid on;
saveas(gcf, 'LOW_IF_CW_1target_plots_IQ/Captured_Signal_zoom.png');


%% Processing
matrix_ch1 = reshape(data(1:M*N), N, M).';

%% Perform MS clutter rejection
% Removing mean from each element in matrix
matrix_ch1_ms = matrix_ch1 - mean(matrix_ch1, 1);

%% Perform FFT
Y_ch1 = fft(matrix_ch1_ms, 4*N, 2);

Y_dB_ch1 = 20*log10(abs(fftshift(Y_ch1,2)));

%% Normalization 1
max_ch1 = max(Y_dB_ch1(:));

% Subtracting the maximum value
Y_dB_ch1_norm1 = Y_dB_ch1 - max_ch1;

%% Plotting

f_axis = linspace(-Fs/2, Fs/2, 4*N);

v_axis = f_axis*c/(fc*2);
t_axis = linspace(1, Tp*M, M);

figure(3)
imagesc(v_axis, t_axis, Y_dB_ch1_norm1, [-55 0])
colorbar;
xlim([-5 5]);
xlabel('Velocity (m/sec)');
ylabel('Time (sec)');
title('Pulse Time Tp = 0.1s, Center Frequency fc = 2.445 GHz Normalization 1');
saveas(gcf, 'LOW_IF_CW_1target_plots_IQ/Norm1.png');

% Velocity and Frequency Resolution
freq_resolution = 1 / Tp;
vel_resolution = c / (2 * fc *Tp);


%% Multi-Target Velocity vs Time with Conditional num_targets
num_targets = 1;  % choose 1 or 2

threshold_dB = -30;  % detection threshold
memory_length = 5;   % number of frames to "remember" a target
max_jump = 2;      % max velocity jump (m/s) allowed for association
MinPeakDistance = 25;

strongest_velocity = zeros(M, num_targets);

% Parameters for tracking
miss_count_single = 0;
miss_count = zeros(M,2);  % counters for missed detections

% When updating tracks, add a "jump condition"
allow_jump = @(new_val, old_val) (old_val == 0) || (abs(new_val - old_val) < max_jump);

for k = 1:M
    % Find all peaks first
    [pks, locs] = findpeaks(Y_dB_ch1_norm1(k,:), 'MinPeakDistance', MinPeakDistance);

    % Keep only peaks above threshold
    valid_idx = pks > threshold_dB;
    pks = pks(valid_idx);
    locs = locs(valid_idx);

    % Convert to velocities
    new_vels = f_axis(locs) * c / (2*fc);

    if num_targets == 1
        if isempty(new_vels)
            if miss_count_single < memory_length
                strongest_velocity(k,1) = strongest_velocity(max(k-1,1),1); % hold previous
                miss_count_single = miss_count_single + 1;
            else
                strongest_velocity(k,1) = 0;
                miss_count_single = 0;
            end
        else
            [~, idx] = max(pks); % strongest peak
            candidate = new_vels(idx);

            if allow_jump(candidate, strongest_velocity(max(k-1,1),1))
                strongest_velocity(k,1) = candidate;
            else
                strongest_velocity(k,1) = strongest_velocity(max(k-1,1),1); % hold previous if jump too big
            end
            miss_count_single = 0; % reset since we got a detection
        
        
        end

    elseif num_targets == 2
        if k == 1
            % Initialization frame: assign directly
            if isempty(new_vels)
                strongest_velocity(k,:) = 0;
            else
                % Take up to 2 strongest peaks
                [~, idx_sort] = sort(pks,'descend');
                top_idx = idx_sort(1:min(2,length(idx_sort)));
                strongest_velocity(k,1:length(top_idx)) = new_vels(top_idx);
            end

        else
            if length(new_vels) >= 1
                % Sort by strength, keep up to num_targets
                [~, idx_sort] = sort(pks, 'descend');
                new_vels = new_vels(idx_sort(1:min(num_targets, length(idx_sort))));
            
                % Build cost matrix (detections x tracks)
                cost_matrix = abs(new_vels(:) - strongest_velocity(k-1,:));
            
                % Solve assignment with cutoff = max_jump
                [assignments, ~] = matchpairs(cost_matrix, max_jump);
            
                % Start with previous values
                strongest_velocity(k,:) = strongest_velocity(k-1,:);
            
                % Update assigned tracks
                for a = 1:size(assignments,1)
                    new_idx = assignments(a,1); % detection index
                    old_idx = assignments(a,2); % track index
                    candidate = new_vels(new_idx);
            
                    % Only update if jump is allowed
                    if abs(candidate - strongest_velocity(k-1, old_idx)) <= max_jump
                        strongest_velocity(k, old_idx) = candidate;
                        miss_count(k, old_idx) = 0;
                    else
                        % Too big jump: keep old value, increase miss count
                        strongest_velocity(k, old_idx) = strongest_velocity(k-1, old_idx);
                        miss_count(k, old_idx) = miss_count(k-1, old_idx) + 1;
                    end
                end
            
                % Handle unassigned tracks (no valid detection nearby)
                assigned_tracks = assignments(:,2);
                unassigned = setdiff(1:num_targets, assigned_tracks);
                for u = unassigned
                    strongest_velocity(k,u) = strongest_velocity(k-1,u);
                    miss_count(k,u) = miss_count(k-1,u) + 1;
                end
            
            else
                % No detections: carry old values and increment miss count
                strongest_velocity(k,:) = strongest_velocity(k-1,:);
                miss_count(k,:) = miss_count(k-1,:) + 1;
            end
            
            % Reset tracks if missed too long
            for ti = 1:num_targets
                if miss_count(k,ti) > memory_length
                    strongest_velocity(k,ti) = 0;
                end
            end


        end

    else
        error('num_targets must be 1 or 2');
    end
end


% Plot velocities
figure(5);
plot(t_axis, strongest_velocity(:,1), '-o', 'DisplayName','Target 1'); hold on;
if num_targets == 2
    plot(t_axis, strongest_velocity(:,2), '-s', 'DisplayName','Target 2');
end
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title(['Velocity vs Time for ', num2str(num_targets), ' Strongest Scatterer(s)']);
legend;
grid on;
saveas(gcf, 'LOW_IF_CW_1target_plots_IQ/Velocity_vs_Time.png');




%% Trend Curve

window = 20;  % adjust (number of frames to average)
vel_smooth1 = movmean(strongest_velocity(:,1), window, 'omitnan');

% order = 5;   % polynomial order
% framelen = 9; % must be odd
% vel_smooth1 = sgolayfilt(strongest_velocity(:,1), order, framelen);

% alpha = 0.5; % 0–1, smaller = smoother
% vel_smooth1 = filter(alpha, [1 alpha-1], strongest_velocity(:,1));

hold on;
plot(t_axis, vel_smooth1, 'b-', 'LineWidth', 2, 'DisplayName','Target 1 (window)');

if num_targets==2
    vel_smooth2 = movmean(strongest_velocity(:,2), window, 'omitnan');
    % vel_smooth2 = sgolayfilt(strongest_velocity(:,2), order, framelen);
    % vel_smooth2 = filter(alpha, [1 alpha-1], strongest_velocity(:,2));
    
    plot(t_axis, vel_smooth2, 'r-', 'LineWidth', 2, 'DisplayName','Target 2 (window)');
end

legend;
saveas(gcf, 'LOW_IF_CW_1target_plots_IQ/Velocity_vs_Time_curves.png');


%% Final Combined Plot: Spectrogram + Velocity Tracks (same axes)

figure(6);
imagesc(v_axis, t_axis, Y_dB_ch1_norm1, [-55 0]);  % spectrogram background
colorbar;
xlim([-5 5]);
xlabel('Velocity (m/sec)');
ylabel('Time (sec)');
title('Spectrogram with Target Velocities');
set(gca,'YDir','reverse');  % make time increase upward
hold on;

% Plot tracked velocities over spectrogram (swap axes to match imagesc)
plot(strongest_velocity(:,1), t_axis, 'b-o', 'LineWidth', 1.5, 'DisplayName','Target 1');
if num_targets == 2
    plot(strongest_velocity(:,2), t_axis, 'r-s', 'LineWidth', 1.5, 'DisplayName','Target 2');
end

% Also plot smoothed velocities
plot(vel_smooth1, t_axis, 'w-', 'LineWidth', 2, 'DisplayName','Target 1 (smoothed)');
if num_targets == 2
    plot(vel_smooth2, t_axis, 'w-', 'LineWidth', 2, 'DisplayName','Target 2 (smoothed)');
end

saveas(gcf, 'LOW_IF_CW_1target_plots_IQ/Spec_with_Velocities.png');
