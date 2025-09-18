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
xlim([0 30]);
xlabel('Velocity (m/sec)');
ylabel('Time (sec)');
title('Pulse Time Tp = 0.1s, Center Frequency fc = 2.445 GHz Normalization 1');
saveas(gcf, 'CW_plots/Norm1.png');

figure(2)
imagesc(v_axis, t_axis, Y_dB_ch1_norm2, [-10 0])
colorbar;
xlim([0 30]);
xlabel('Velocity (m/sec)');
ylabel('Time (sec)');
title('Pulse Time Tp = 0.1s, Center Frequency fc = 2.445 GHz Normalization 2');
saveas(gcf, 'CW_plots/Norm2.png');

% Velocity and Frequency Resolution
% freq_resolution = 1 / Tp;
% vel_resolution = c / (2 * fc *Tp);


