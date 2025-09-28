%% Velocity Measurements - CW Mode - Single and Multiple Targets
% 
% Group 7
% Review Meeting 3
% 23/09/2025
% 

%% Cleaning
clear;
clc;
close all;

%% Read data
[data, Fs] = audioread("Test_files/cw_18_t2.wav");
data=data(:,2);

Tp = 0.1;
fc = 2.445e9;
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
saveas(gcf, 'CW_plots/Captured_Signal_total.png');

figure(2);
plot(t, data);
xlabel('Time (s)');
ylabel('Amplitude');
title('Captured Radar Signal');
xlim([5 5.05]);
grid on;
saveas(gcf, 'CW_plots/Captured_Signal_zoom.png');


%% Processing
% matrix_ch1 = zeros(M,N);
% sum_ch1 = 0;

% for i = 0:M-1
%     for j = 1:N
%         if i*N+j <= length(data)
%             matrix_ch1(i+1,j) = data(i*N+j,1); % Inserting element in matrix
%             sum_ch1 = sum_ch1 + matrix_ch1(i+1,j); % Calculating partial sum to get mean
%         end
%     end
% end

matrix_ch1 = reshape(data(1:M*N), N, M).';

%% Perform MS clutter rejection

% Removing mean from each element in matrix
% matrix_ch1_ms = matrix_ch1 - mean_ch1;
matrix_ch1_ms = matrix_ch1 - mean(matrix_ch1, 1);

%% Perform FFT
Y_ch1 = fft(matrix_ch1_ms, 4*N, 2);

% Decibel values
Y_dB_ch1 = 20*log10(abs(Y_ch1));

% Calculating Fmax value
Fmax = Fs/2;

% Half of the matrix
Y_dB_ch1_half = Y_dB_ch1(:,1:4*N/2);


%% Normalization 1
max_ch1 = max(Y_dB_ch1_half(:));

% Subtracting the maximum value
Y_dB_ch1_norm1 = Y_dB_ch1_half - max_ch1;


%% Normalization 2

Y_dB_ch1_norm2 = Y_dB_ch1_half - max(Y_dB_ch1_half, [], 2);


%% Plotting

f_axis = linspace(0, Fmax, 4*N/2);
v_axis = f_axis*c/(fc*2);
t_axis = linspace(1, Tp*M, M);

figure(3)
imagesc(v_axis, t_axis, Y_dB_ch1_norm1, [-55 0])
colorbar;
% xlim([0 30]);
xlim([0 10]);
xlabel('Velocity (m/sec)');
ylabel('Time (sec)');
title('Pulse Time Tp = 0.1s, Center Frequency fc = 2.445 GHz Normalization 1');
saveas(gcf, 'CW_plots/Norm1.png');

figure(4)
imagesc(v_axis, t_axis, Y_dB_ch1_norm2, [-10 0])
colorbar;
% xlim([0 30]);
xlim([0 10]);
xlabel('Velocity (m/sec)');
ylabel('Time (sec)');
title('Pulse Time Tp = 0.1s, Center Frequency fc = 2.445 GHz Normalization 2');
saveas(gcf, 'CW_plots/Norm2.png');

% Velocity and Frequency Resolution
freq_resolution = 1 / Tp;
vel_resolution = c / (2 * fc *Tp);


%% Multi-Target Velocity vs Time with Conditional num_targets
num_targets = 1;  % choose 1 or 2

strongest_velocity = nan(M, num_targets);

threshold_dB = -20;  % detection threshold

% Parameters for tracking
miss_count_single = 0;
memory_length = 3;   % number of frames to "remember" a target
max_jump = 2;        % max velocity jump (m/s) allowed for association
miss_count = zeros(M,2);  % counters for missed detections
% When updating tracks, add a "jump condition"
allow_jump = @(new_val, old_val) (old_val == 0) || (abs(new_val - old_val) < max_jump);

for k = 1:M
    % Find all peaks first
    [pks, locs] = findpeaks(Y_dB_ch1_norm1(k,:), 'MinPeakDistance', 25);

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
            strongest_velocity(k,1) = new_vels(idx);
        end

    elseif num_targets == 2
        if k == 1
            % Initialization frame → assign directly
            if isempty(new_vels)
                strongest_velocity(k,:) = 0;
            else
                % Take up to 2 strongest peaks
                [~, idx_sort] = sort(pks,'descend');
                top_idx = idx_sort(1:min(2,length(idx_sort)));
                strongest_velocity(k,1:length(top_idx)) = new_vels(top_idx);
            end

        else
            if length(new_vels) >= 2
                % Take 2 strongest peaks
                [~, idx_sort] = sort(pks,'descend');
                new_vels = new_vels(idx_sort(1:2));

                % Compute assignment costs
                dist11 = abs(new_vels(1) - strongest_velocity(k-1,1));
                dist12 = abs(new_vels(1) - strongest_velocity(k-1,2));
                dist21 = abs(new_vels(2) - strongest_velocity(k-1,1));
                dist22 = abs(new_vels(2) - strongest_velocity(k-1,2));

                cost_assign1 = dist11 + dist22; % new1→old1, new2→old2
                cost_assign2 = dist12 + dist21; % new1→old2, new2→old1

                if cost_assign1 <= cost_assign2
                    if allow_jump(new_vels(1), strongest_velocity(k-1,1))
                        strongest_velocity(k,1) = new_vels(1);
                    else
                        strongest_velocity(k,1) = strongest_velocity(k-1,1); % hold old
                    end
                    if allow_jump(new_vels(2), strongest_velocity(k-1,2))
                        strongest_velocity(k,2) = new_vels(2);
                    else
                        strongest_velocity(k,2) = strongest_velocity(k-1,2);
                    end
                else
                    if allow_jump(new_vels(2), strongest_velocity(k-1,1))
                        strongest_velocity(k,1) = new_vels(2);
                    else
                        strongest_velocity(k,1) = strongest_velocity(k-1,1);
                    end
                    if allow_jump(new_vels(1), strongest_velocity(k-1,2))
                        strongest_velocity(k,2) = new_vels(1);
                    else
                        strongest_velocity(k,2) = strongest_velocity(k-1,2);
                    end
                end

                miss_count(k,:) = 0; % reset counters

            elseif isscalar(new_vels)
                % One detection → assign to closest previous track
                d1 = abs(new_vels - strongest_velocity(k-1,1));
                d2 = abs(new_vels - strongest_velocity(k-1,2));

                if d1 < d2 && d1 < max_jump
                    strongest_velocity(k,1) = new_vels;
                    strongest_velocity(k,2) = strongest_velocity(k-1,2); % keep old
                    miss_count(k,2) = miss_count(max(k-1,1),2)+1;
                elseif d2 < max_jump
                    strongest_velocity(k,2) = new_vels;
                    strongest_velocity(k,1) = strongest_velocity(k-1,1);
                    miss_count(k,1) = miss_count(max(k-1,1),1)+1;
                else
                    % Detection doesn't match → carry old velocities
                    strongest_velocity(k,:) = strongest_velocity(k-1,:);
                    miss_count(k,:) = miss_count(max(k-1,1),:)+1;
                end

            else
                % No detections → hold previous velocities (if no detection<3 points)
                strongest_velocity(k,:) = strongest_velocity(k-1,:);
                miss_count(k,:) = miss_count(max(k-1,1),:)+1;
            end

            % Reset to 0 if missing too long
            for ti = 1:2
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
saveas(gcf, 'CW_plots/Velocity_vs_Time.png');




%% Trend Curve

window = 5;  % adjust (number of frames to average)
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
saveas(gcf, 'CW_plots/Velocity_vs_Time_curves.png');


%% Final Combined Plot: Spectrogram + Velocity Tracks (same axes)

figure(6);
imagesc(v_axis, t_axis, Y_dB_ch1_norm1, [-55 0]);  % spectrogram background
colorbar;
xlim([0 10]);
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
    plot(vel_smooth2, t_axis, 'm-', 'LineWidth', 2, 'DisplayName','Target 2 (smoothed)');
end

saveas(gcf, 'CW_plots/Spec_with_Velocities.png');
