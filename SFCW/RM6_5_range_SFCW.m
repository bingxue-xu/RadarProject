%% Range Measurement – SFCW (Two-Step Mode)
clear; clc; close all;

%% ---- USER PARAMETERS ----
c  = 3e8;              % speed of light [m/s]
f1 = 5.80e9;           % first frequency [Hz]
df = 50e6;             % step size (Hz)
f2 = f1 + df;          % second frequency
N  = 2;                % number of frequency steps
% Fs = 2e6;              % sample rate (adjust to your SDR)
Ts = 0.001;             % dwell time per frequency step (sec)
% M  = 200;              % number of 2-step cycles

%% ---- LOAD DATA ----
[data_I, Fs] = audioread('cw_im_2steps/cw_re_2steps.wav');
[data_Q, ~]  = audioread('cw_im_2steps/cw_im_2steps.wav');

rx = data_I + 1j*data_Q;

% (Optional) known square wave sync channel if recorded
% sync = data(:,3) > 0;  % threshold detect if exists

t = (0:length(rx)-1)/Fs; 

%% ---- SEGMENT DATA BY FREQUENCY STEP ----
samples_per_step = round(Ts*Fs);
samples_per_cycle = N * samples_per_step;

M = floor(length(rx)/N/2);

S1 = zeros(M,1);
S2 = zeros(M,1);

for m = 1:M
    base = (m-1)*samples_per_cycle + 1;
    idx1 = base : base + samples_per_step - 1;
    idx2 = base + samples_per_step : base + 2*samples_per_step - 1;

    if max(idx2) > length(rx)
        break
    end

    % Average IQ within each dwell to get one complex value per frequency
    S1(m) = mean(rx(idx1));
    S2(m) = mean(rx(idx2));
end

S1 = S1 - mean(S1);
S2 = S2 - mean(S2);   

%% ---- PHASE DIFFERENCE & RANGE ESTIMATION ----
dphi = angle(S2 .* conj(S1))+pi;          % (-0, 2*pi]
R = (c * dphi) / (4 * pi * df);        % meters, within ±1.5 m
Rmax = c / (4 * df);                   % 1.5 m for df=50 MHz

fprintf('Theoretical unambiguous range = %.2f m\n', Rmax);

%% ---- TIME AXIS ----
time_axis = (0:length(R)-1) * (2*Ts);  % one cycle = 2*Ts

%% ---- PLOT RANGE OVER TIME ----
figure();
subplot(1,2,1)
plot(time_axis, R, 'b.-','LineWidth',1.2);
xlabel('Time (s)');
xlim([0,20]);
ylabel('Estimated Range (m)');
title('SFCW Two-Step Range vs. Time (within 1.5 m window)');
grid on;

%% ---- PLOT RANGE OVER TIME ----
subplot(1,2,2)
plot(time_axis, R, 'b.-','LineWidth',1.2);
xlabel('Time (s)');
xlim([10,11]);
ylabel('Estimated Range (m)');
title('SFCW Two-Step Range vs. Time (within 1.5 m window)');
grid on;

% %% ---- OPTIONAL: unwrap if small motion causes ±π jumps ----
% R_unwrapped = unwrap(dphi) * (c / (4*pi*df));
% figure;
% plot(time_axis, R_unwrapped, 'r.-');
% xlabel('Time (s)');
% ylabel('Unwrapped Range (m)');
% title('Unwrapped Range (if motion crosses ±π)');
% grid on;

figure;
plot(time_axis, mod(dphi, 2*pi), 'b.-');
xlabel('Time (s)');
xlim([0,20]);
ylabel('Phase Difference (rad)');
title('Phase Difference (wrapped 0–2π)');
grid on;
ylim([0 2*pi]);

%% ---- Build IQ-Time Matrix ----
IQ_matrix = [S1, S2];   % Each column = one frequency step (f1, f2)
M = length(S1);

%% ---- Plot Correctly ----
figure;

subplot(2,1,1);
plot(real(S1), 'b'); hold on;
plot(real(S2), 'r');
title('Before Clutter Rejection');
xlabel('Cycle index (time)');
xlim([0, 10000]);
ylabel('Amplitude');
legend('f1','f2');
grid on;

subplot(2,1,2);
imagesc(1:2, 1:M, abs(IQ_matrix));    % correct x-axis = 1..2
xlabel('Frequency Step (f1 → f2)');
ylabel('Cycle Index (Time)');
title('IQ-Time Magnitude Matrix');
colorbar;


T = table((1:10)', real(S1(1:10)), imag(S1(1:10)), ...
           real(S2(1:10)), imag(S2(1:10)), ...
           'VariableNames', {'Cycle','I_f1','Q_f1','I_f2','Q_f2'});
disp(T);

dphi = angle(S2 .* conj(S1));
R = (c * dphi) / (4 * pi * df);

T_range = table((1:M)', R, dphi, ...
    'VariableNames', {'Cycle','Range_m','Phase_rad'});
disp(T_range(1:10,:));
