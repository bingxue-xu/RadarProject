%% Velocity Measurements - CW Doppler Radar
% 
% Velocity Algorithm development for CW Doppler Radar.
% Effect of variation in Tp (above and below the nominal value) on the 
% velocity resolution and velocity accuracy.
% Variation in results with and without MS clutter rejection.
% Variation in results utilizing more or less zero padding for FFT.
% Variation in results with the two normalization methods.
% 
% Group 7
% 10/09/2025

%% Cleaning
clear;
clc;

%% Read data
[data, Fs] = audioread("Velocity_Test_File2.m4a");

Tp = 0.1; % For subtask: above and below the nominal value
fc = 2.445e9;
c = 3e8;
N = Tp*Fs;
M = ceil(length(data)/N);

matrix_ch1 = zeros(M,N);
matrix_ch2 = zeros(M,N);
sum_ch1 = 0;
sum_ch2 = 0;

for i = 0:M-1
    for j = 1:N
        if i*N+j <= length(data)
            matrix_ch1(i+1,j) = data(i*N+j,1); % Inserting element in matrix
            sum_ch1 = sum_ch1 + matrix_ch1(i+1,j); % Calculating partial sum to get mean

            matrix_ch2(i+1,j) = data(i*N+j,2);
            sum_ch2 = sum_ch2 + matrix_ch2(i+1,j);
        end
    end
end

% Evaluate mean
mean_ch1 = sum_ch1 / length(data);
mean_ch2 = sum_ch2 / length(data);

%% Perform MS clutter rejection

% Removing mean from each element in matrix
matrix_ch1_ms = matrix_ch1 - mean_ch1;
matrix_ch2_ms = matrix_ch2 - mean_ch2;

%% Perform FFT
Y_ch1_raw = fft(matrix_ch1, 4*N, 2);
Y_ch2_raw = fft(matrix_ch2, 4*N, 2);

Y_ch1 = fft(matrix_ch1_ms, 4*N, 2);
Y_ch2 = fft(matrix_ch2_ms, 4*N, 2);

Y_ch1_8zp = fft(matrix_ch1_ms, 8*N, 2);
Y_ch1_0zp = fft(matrix_ch1_ms, [], 2);

% Decibel values
Y_dB_ch1_raw = 20*log10(abs(Y_ch1_raw));
Y_dB_ch2_raw = 20*log10(abs(Y_ch2_raw));

Y_dB_ch1 = 20*log10(abs(Y_ch1));
Y_dB_ch2 = 20*log10(abs(Y_ch2));

Y_dB_ch1_8zp = 20*log10(abs(Y_ch1_8zp));
Y_dB_ch1_0zp = 20*log10(abs(Y_ch1_0zp));

% Calculating Fmax value
Fmax = Fs/2;

% Half of the matrix
Y_dB_ch1_half_raw = Y_dB_ch1_raw(:,1:2*N);
Y_dB_ch2_half_raw = Y_dB_ch2_raw(:,1:2*N);

Y_dB_ch1_half = Y_dB_ch1(:,1:2*N);
Y_dB_ch2_half = Y_dB_ch2(:,1:2*N);

Y_dB_ch1_half_8zp = Y_dB_ch1_8zp(:,1:4*N);
Y_dB_ch1_half_0zp = Y_dB_ch1_0zp(:,1:N/2);


%% Normalization 1

max_ch1_raw = max(Y_dB_ch1_half_raw(:));
max_ch2_raw = max(Y_dB_ch2_half_raw(:));

max_ch1 = max(Y_dB_ch1_half(:));
max_ch2 = max(Y_dB_ch2_half(:));

max_ch1_8zp = max(Y_dB_ch1_half_8zp(:));
max_ch1_0zp = max(Y_dB_ch1_half_0zp(:));

% Subtracting the maximum value
Y_dB_ch1_norm1_raw = Y_dB_ch1_half_raw - max_ch1_raw;
Y_dB_ch2_norm1_raw = Y_dB_ch2_half_raw - max_ch2_raw;

Y_dB_ch1_norm1 = Y_dB_ch1_half - max_ch1;
Y_dB_ch2_norm1 = Y_dB_ch2_half - max_ch2;

Y_dB_ch1_norm1_8zp = Y_dB_ch1_half_8zp - max_ch1_8zp;
Y_dB_ch1_norm1_0zp = Y_dB_ch1_half_0zp - max_ch1_0zp;

%% Normalization 2
Y_dB_ch1_norm2_raw = zeros(M,2*N);
Y_dB_ch2_norm2_raw = zeros(M,2*N);

Y_dB_ch1_norm2 = zeros(M,2*N);
Y_dB_ch2_norm2 = zeros(M,2*N);

% Find the maximum value in each row
max_row_ch1 = max(Y_dB_ch1_half, [], 2);
max_row_ch1_raw = max(Y_dB_ch1_half_raw, [], 2);

% Subtract that value in each element of the row
for i = 1:length(max_row_ch1)
    Y_dB_ch1_norm2(i,:) = Y_dB_ch1_half(i,:) - max_row_ch1(i);
    Y_dB_ch1_norm2_raw(i,:) = Y_dB_ch1_half_raw(i,:) - max_row_ch1_raw(i);
end

max_row_ch2 = max(Y_dB_ch2_half, [], 2);

for i = 1:length(max_row_ch2)
    Y_dB_ch2_norm2(i,:) = Y_dB_ch2_half(i,:) - max_row_ch2(i);
end


%% Plotting

f_axis = linspace(0, Fmax, 2*N);
v_axis = f_axis*c/(fc*2);
t_axis = linspace(1, Tp*M, M);

figure(1)
imagesc(v_axis, t_axis, Y_dB_ch1_norm1, [-55 0])
colorbar;
xlim([0 30]);
xlabel('Velocity (m/sec)');
ylabel('Time (sec)');
title('Pulse Time Tp = 0.1s, Center Frequency fc = 2.445 GHz Normalization 1');
% saveas(gcf, 'Norm1_Tp_0_1.png');

figure(2)
imagesc(v_axis, t_axis, Y_dB_ch1_norm2, [-10 0])
colorbar;
xlim([0 30]);
xlabel('Velocity (m/sec)');
ylabel('Time (sec)');
title('Pulse Time Tp = 0.1s, Center Frequency fc = 2.445 GHz Normalization 2');
% saveas(gcf, 'Norm2_Tp_0_1.png');


% Without MS clutter rejection

% Raw data
figure(3)
imagesc(v_axis, t_axis, Y_dB_ch1_norm1_raw, [-55 0]);
colorbar;
xlim([0 30]);
xlabel('Velocity (m/sec)');
ylabel('Time (sec)');
title('Pulse Time Tp = 0.1s, without MS, Center Frequency fc = 2.445 GHz Normalization 1');
% saveas(gcf, 'Norm1_Tp_0_1_no_MS_mscomp.png');

% MS
figure(4)
imagesc(v_axis, t_axis, Y_dB_ch1_norm1, [-55 0])
colorbar;
xlim([0 30]);
xlabel('Velocity (m/sec)');
ylabel('Time (sec)');
title('Pulse Time Tp = 0.1s, with MS, Center Frequency fc = 2.445 GHz Normalization 1');
% saveas(gcf, 'No_Norm_Tp_0_1_MS_mscomp.png');


% More and less zero padding

figure(5)
imagesc(v_axis, t_axis, Y_dB_ch1_norm1_8zp, [-55 0])
colorbar;
xlim([0 30]);
xlabel('Velocity (m/sec)');
ylabel('Time (sec)');
title('Pulse Time Tp = 0.1s, with MS and zp 8N, Center Frequency fc = 2.445 GHz Normalization 1');
% saveas(gcf, 'No_Norm_8zp.png');

figure(6)
imagesc(v_axis, t_axis, Y_dB_ch1_norm1_0zp, [-55 0])
colorbar;
xlim([0 30]);
xlabel('Velocity (m/sec)');
ylabel('Time (sec)');
title('Pulse Time Tp = 0.1s, with MS and no zp, Center Frequency fc = 2.445 GHz Normalization 1');
% saveas(gcf, 'No_Norm_0zp.png');

% Raw data
figure(7)
imagesc(v_axis, t_axis, Y_dB_ch1_norm1_raw, [-55 0])
colorbar;
xlim([0 30]);
xlabel('Velocity (m/sec)');
ylabel('Time (sec)');
title('Pulse Time Tp = 0.1s, raw data, Center Frequency fc = 2.445 GHz Normalization 1');
% saveas(gcf, 'Raw_data_norm1.png');

% Velocity and Frequency Resolution
freq_resolution = 1 / Tp;
vel_resolution = c / (2 * fc *Tp);


