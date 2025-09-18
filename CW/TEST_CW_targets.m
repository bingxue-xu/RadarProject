%% Velocity Measurements - CW Doppler Radar

%% Cleaning
clear;
clc;

%% Read data
[data, Fs] = audioread("Velocity_Test_CW_1.wav");
data=data(:,2);

Tp = 0.1;
fc = 2.445e9;
c = 3e8;
N = Tp*Fs;

M = floor(length(data)/N);

%% Plot captured data
t = (0:length(data)-1)/Fs;  % time axis in seconds

figure(3);
plot(t, data);
xlabel('Time (s)');
ylabel('Amplitude');
title('Captured Radar Signal');
grid on;
saveas(gcf, 'CW_plots/Captured_Signal.png');

figure(4);
plot(t, data);
xlabel('Time (s)');
ylabel('Amplitude');
title('Captured Radar Signal');
xlim([5.9 6.1]);
grid on;
saveas(gcf, 'CW_plots/Captured_Signal.png');


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

figure(1)
imagesc(v_axis, t_axis, Y_dB_ch1_norm1, [-55 0])
colorbar;
% xlim([0 30]);
xlim([0 10]);
xlabel('Velocity (m/sec)');
ylabel('Time (sec)');
title('Pulse Time Tp = 0.1s, Center Frequency fc = 2.445 GHz Normalization 1');
saveas(gcf, 'CW_plots/Norm1.png');

figure(2)
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
num_targets = 2;  % choose 1 or 2

strongest_velocity = nan(M, num_targets);

threshold_dB = -20;  % detection threshold

for k = 1:M
    % Find all peaks first
    [pks, locs] = findpeaks(Y_dB_ch1_norm1(k,:), 'MinPeakDistance', 25);
    
    % Keep only peaks above threshold
    valid_idx = pks > threshold_dB;
    pks = pks(valid_idx);
    locs = locs(valid_idx);
    
    if isempty(pks)
        % No target detected → assign velocity = 0
        strongest_velocity(k,:) = 0;
        continue;
    end
    
    % Sort peaks by amplitude
    [~, idx_sort] = sort(pks, 'descend');
    
    if num_targets == 1
        % Take only the strongest peak
        top_idx = locs(idx_sort(1));
        strongest_velocity(k,1) = f_axis(top_idx) * c / (2*fc);
        
    elseif num_targets == 2
        % Take up to two strongest peaks
        top_idx = locs(idx_sort(1:min(2,length(locs))));
        strongest_velocity(k,1:length(top_idx)) = f_axis(top_idx) * c / (2*fc);
        
        % If only one peak was valid, fill second with 0
        if length(top_idx) < 2
            strongest_velocity(k,2) = 0;
        end
        
    else
        error('num_targets must be 1 or 2');
    end
end

% Plot velocities
figure(3);
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




