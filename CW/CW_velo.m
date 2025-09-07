%% Cleaning
clear;
clc;

%% Get datas
[data, Fs] = audioread("Velocity_Test_File2.wav");

Tp = 0.1;
fc = 2.445e9;
c = 3e8;
N = Tp*Fs;
M = ceil(length(data)/N);

matrix_ch1 = zeros(M, N); % Preallocate the matrix for efficiency
matrix_ch2 = zeros(M,N);
sum_ch1 = 0;
sum_ch2 = 0;

for i = 0:M-1
    for j = 1:N
        if i*N+j <= length(data)
            % Inserting element in matrix
            matrix_ch1(i+1,j) = data(i*N+j,1);
            % Calculatin partial sum to get mean
            sum_ch1 = sum_ch1 + matrix_ch1(i+1,j);
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
matrix_ch1 = matrix_ch1 - mean_ch1;
matrix_ch2 = matrix_ch2 - mean_ch2;

%% Perform FFT
Y_ch1 = fft(matrix_ch1, 4*N, 2);
Y_ch2 = fft(matrix_ch2, 4*N, 2);

% Decibel values
Y_dB_ch1 = 20*log10(abs(Y_ch1));
Y_dB_ch2 = 20*log10(abs(Y_ch2));

% Calculating Fmax value
Fmax = Fs/2;
Y_dB_ch1_half = Y_dB_ch1(:,1:2*N);
Y_dB_ch2_half = Y_dB_ch2(:,1:2*N);

%% Normalization 1

% Find the maximum for each row, and then the maximum among these values
max_ch1 = max(Y_dB_ch1_half(:));
max_ch2 = max(Y_dB_ch2_half(:));
% Subtracting the maximum value
Y_dB_ch1_norm1 = Y_dB_ch1_half - max_ch1;
Y_dB_ch2_norm1 = Y_dB_ch2_half - max_ch2;

%% Normalization 2
Y_dB_ch1_norm2 = zeros(M,2*N);
Y_dB_ch2_norm2 = zeros(M,2*N);
% Find the maximum value in each row
max_row_ch1 = max(Y_dB_ch1_half, [], 2);
% Subtract that value in each element of the row
for i = 1:length(max_row_ch1)
    Y_dB_ch1_norm2(i,:) = Y_dB_ch1_half(i,:) - max_row_ch1(i);
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

figure(2)
imagesc(v_axis, t_axis, Y_dB_ch1_norm2, [-10 0])
colorbar;
xlim([0 30]);
xlabel('Velocity (m/sec)');
ylabel('Time (sec)');
title('Pulse Time Tp = 0.1s, Center Frequency fc = 2.445 GHz Normalization 2');