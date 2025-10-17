%% Range Measurement – SFCW (Two-Step Mode)
clear; clc; close all;

%% ---- USER PARAMETERS ----
c  = 3e8;              % speed of light [m/s]
f1 = 5.80e9;           % first frequency [Hz]
% df = 50e6;             % step size (Hz)
df = 25e6;
f2 = f1 + df;          % second frequency
f3 = f1 + 2*df;
f4 = f1 + 3*df;

N  = 4;                % number of frequency steps
% Fs = 2e6;              % sample rate (adjust to your SDR)
% Ts = 0.05;             % dwell time per frequency step (sec)
Ts = 1e-3;
M  = 5000;              % number of 2-step cycles

%% ---- LOAD DATA ----
[data_I, Fs] = audioread('cw_im_4steps_lasttest/cw_re_4steps_lasttest_3m.wav');
[data_Q, ~]  = audioread('cw_im_4steps_lasttest/cw_im_4steps_lasttest_3m.wav');

rx = data_I + 1j*data_Q;

% (Optional) known square wave sync channel if recorded
% sync = data(:,3) > 0;  % threshold detect if exists

t = (0:length(rx)-1)/Fs; 

%% ---- SEGMENT DATA BY FREQUENCY STEP ----
samples_per_step = round(Ts*Fs);
samples_per_cycle = N * samples_per_step;

S1 = zeros(M,1);
S2 = zeros(M,1);
S3 = zeros(M,1);
S4 = zeros(M,1);

for m = 1:M
    base = (m-1)*samples_per_cycle + 1;
    idx1 = base : base + samples_per_step - 1;
    idx2 = base + samples_per_step : base + 2*samples_per_step - 1;

    idx3 = base + 2*samples_per_step : base + 3*samples_per_step - 1;
    idx4 = base + 3*samples_per_step : base + 4*samples_per_step - 1;

    if max(idx4) > length(rx)
        break
    end

    % Average IQ within each dwell to get one complex value per frequency
    S1(m) = mean(rx(idx1));
    S2(m) = mean(rx(idx2));
    S3(m) = mean(rx(idx3));
    S4(m) = mean(rx(idx4));
end

%% ---- PHASE DIFFERENCE & RANGE ESTIMATION ----
dphi_1 = angle(S2 .* conj(S1));          % (-pi, pi]
dphi_2 = angle(S3 .* conj(S2));          % (-pi, pi]
dphi_3 = angle(S4 .* conj(S3));          % (-pi, pi]

% dphi = [dphi_1 dphi_2 dphi_3];
dphi = reshape([dphi_1(:) dphi_2(:) dphi_3(:)].', [], 1);

%

R_1 = (c * dphi_1) / (4 * pi * df);        % meters, within ±1.5 m
R_2 = (c * dphi_2) / (4 * pi * df);        % meters, within ±1.5 m
R_3 = (c * dphi_3) / (4 * pi * df);        % meters, within ±1.5 m

% R = [R_1 R_2 R_3];
R = reshape([R_1(:) R_2(:) R_3(:)].', [], 1);

% Rmax = c / (4 * df);                   % 1.5 m for df=50 MHz
Rmax = c / (2 * df);                   % 1.5 m for df=50 MHz

fprintf('Theoretical unambiguous range = %.2f m\n', Rmax);

%% ---- TIME AXIS ----
% time_axis = (0:length(R)-1) * (N*Ts);  % one cycle = 2*Ts
time_axis = (0:length(R)-1) * Ts;


%% ---- PLOT RANGE OVER TIME ----
figure;
plot(time_axis, R, 'b.-','LineWidth',1.2);
xlabel('Time (s)');
ylabel('Estimated Range (m)');
title('SFCW Two-Step Range vs. Time (within 1.5 m window)');
grid on;

%% ---- OPTIONAL: unwrap if small motion causes ±π jumps ----
R_unwrapped = unwrap(dphi) * (c / (4*pi*N*df));
figure;
plot(time_axis, R_unwrapped, 'r.-');
xlabel('Time (s)');
ylabel('Unwrapped Range (m)');
title('Unwrapped Range (if motion crosses ±π)');
grid on;

%% Phase plot
figure;
plot(time_axis, dphi, 'r.-');
xlabel('Time (s)');
ylabel('Phase (rad)');
title('Unwrapped Range (if motion crosses ±π)');
grid on;
figure;
plot(time_axis, abs(unwrap(dphi)), 'r.-');
xlabel('Time (s)');
ylabel('Unwrapped Phase (rad)');
title('Unwrapped Range (if motion crosses ±π)');
grid on;

