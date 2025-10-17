%% Range Measurement – SFCW (Two-Step Mode)
clear; clc; close all;

%% ---- USER PARAMETERS ----
c  = 3e8;              % speed of light [m/s]
f1 = 5.80e9;           % first frequency [Hz]
% df = 50e6;             % step size (Hz)
df = 50e6;
f2 = f1 + df;          % second frequency
f3 = f1 + 2*df;
f4 = f1 + 3*df;

N  = 2;                % number of frequency steps
% Fs = 2e6;              % sample rate (adjust to your SDR)
% Ts = 0.05;             % dwell time per frequency step (sec)
Ts = 100e-3;
% M  = 5000;              % number of 2-step cycles

%% ---- LOAD DATA ----
[data_I, Fs] = audioread('cw_im_2steps/cw_re_2steps.wav');
[data_Q, ~]  = audioread('cw_im_2steps/cw_im_2steps.wav');

rx = data_I + 1j*data_Q;


%% ---- SEGMENT DATA BY FREQUENCY STEP ----
samples_per_step  = round(Ts*Fs);
samples_per_cycle = N * samples_per_step;
M = floor(length(rx)/samples_per_cycle);
L = M*N*samples_per_step;
rx = rx(1:L);

%% matrix
% rx3 = reshape(rx, samples_per_step, N, M);         % [samples, freq, cycles]
% S_matrix = squeeze(mean(rx3,1)).';      % M×N

S_matrix = zeros(M, N);  % Initialize the S_matrix to store averages
for m = 1:M
    base = (m-1)*samples_per_cycle;
    for n = 1: N
        seg = base + (n-1)*samples_per_step + (1: samples_per_step);
        S_matrix(m, n) = mean(rx(seg));  % Store the average for each segment
    end
end

S_matrix = S_matrix - mean(S_matrix);

%% ---- FFT ----
Nfft = 4*N;
R_spec = fft(S_matrix, Nfft, 2);
half = 1:(Nfft/2);
R_mag_half = abs(R_spec(:,half));                   % linear for detection
R_dB_half  = 20*log10(R_mag_half + eps);            % dB for plot

%% two-step plot
dphi12 = angle(S_matrix(:,2).*conj(S_matrix(:,1)));
R12 = (c*dphi12)/(4*pi*df);
Rmax = c / (4 * df);                   % 1.5 m for df=50 MHz
fprintf('Theoretical unambiguous range = %.2f m\n', Rmax);

% map to FFT bin index and overlay
k_est = round(R12 * (2*Nfft*df)/c);                 % 0..Nfft/2-1
k_est = min(max(k_est,0), Nfft/2-1) + 1;            % clamp to valid
R_axis = (half-1) * (c/(2*Nfft*df));
time_axis = (0:size(S_matrix,1)-1) * (N*Ts);

figure; imagesc(R_axis, time_axis, R_dB_half);
axis xy; colorbar; xlabel('Range (m)'); ylabel('Time (s)');
title('SFCW Range–Time (FFT across frequency)'); hold on;
plot(R_axis(k_est), time_axis, 'r.', 'MarkerSize', 10);    % should sit on the ridge
set(gca,'YDir','reverse');



