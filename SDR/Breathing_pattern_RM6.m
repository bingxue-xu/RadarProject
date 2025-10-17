%% Breathing Pattern - CW and Low-IF CW Mode
% 
% Group 7
% Review Meeting 6
% 14/10/2025
% 

clear;
clc;
close all;

%% Parameters
c = 3e8;              % speed of light [m/s]
fc = 5.8e9;           % radar carrier frequency [Hz]
lambda = c / fc;      % wavelength [m]

%% Load Data
[i_data, Fs] = audioread('breathing_samples/low_if_cw_re_breathing_last6.wav');
[q_data, ~]  = audioread('breathing_samples/low_if_cw_im_breathing_last6.wav');

% used to remove some of the initial and final parts
i_data = i_data(1*Fs:(end-Fs*10));
q_data = q_data(1*Fs:(end-Fs*10));

% Use the same number of samples: to check if they have the same
% number of samples to be safer
N = min(length(i_data), length(q_data));
i_data = i_data(1:N);
q_data = q_data(1:N);

analytic_signal = i_data + 1j * q_data;
phase_signal = unwrap(angle(analytic_signal));
t = (0:N-1) / Fs; % time vector

figure(1);
plot(t, phase_signal, 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Phase [rad]');
title('Chest Expansion / Contraction');
grid on;


%% Range Variation (breathing)
range_variation = (c / (4*pi*fc)) * (phase_signal); % meters

figure(2);
plot(t, range_variation*1000, 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Range Variation [mm]');
title('Chest Expansion / Contraction');
grid on;

%% FFT Analysis (Breathing Rate)

L = length(range_variation)*4;
Y = fft(range_variation,L);

Y = Y(1:L/2);
f_axis = linspace(0, Fs/2, L/2).';

% Y = Y(1:L/4);
% f_axis = linspace(0, Fs/2, L/4).';

P = abs(Y);

figure(3);
plot(f_axis, P, 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of Range Variation');
xlim([0 2]);
grid on;

figure(4);
plot(f_axis, 20*log10(P), 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude [dB]');
title('Frequency Spectrum of Range Variation');
xlim([0 2]);
% ylim([-10 60]);
grid on;

P = abs(Y/max(Y));

figure(41);
plot(f_axis, 20*log10(P), 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude Normalized [dB]');
title('Frequency Spectrum of Range Variation');
xlim([0 2]);
% ylim([-60 10]);
grid on;

%% Identify Peaks (Breathing Rate)
[pks, locs] = findpeaks(P, f_axis, 'SortStr', 'descend', 'NPeaks', 2);

disp('Dominant Motion Frequencies (Hz):');
disp(locs);

disp('Breath Rates (Breaths per Minute):');
disp(locs * 60);

% Note: most of the time the 2nd peak shows the real Breathing Rate

