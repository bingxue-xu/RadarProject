%% Range Measurements - SFCW Mode
% 
% Group 7
% Review Meeting 6
% 14/10/2025

clc; close all; clearvars;

%% Settings
% For simulated data: Nf, Ts, df, M
% For real data: Nf, Ts, df

Nf = 2;                  % number of frequency steps per sweep
Ts = 200e-3;             % time of each frequency step
df = 50e6;               % frequency step

M  = 500;                % number of sweeps

f0 = 5.8e9;              % start frequency
c  = 3e8;

%% Real data
[data_I, Fs] = audioread('samples_12_10_indoors/cw_re_200ms_N2_in.wav');
[data_Q, ~]  = audioread('samples_12_10_indoors/cw_im_200ms_N2_in.wav');

outdir = sprintf('Results_temp/Ts%.0fms_Nf%d_df%.0fMHz', Ts*1e3, Nf, df/1e6);
if ~exist(outdir, 'dir')
    mkdir(outdir); 
end

% Pre-processing to avoid spikes
thresh = 0.02;
clean_signal = @(x) interp1(find(abs(x) <= thresh), x(abs(x) <= thresh), 1:length(x), 'linear', 'extrap');
data_I_clean = clean_signal(data_I); data_Q_clean = clean_signal(data_Q);
data_I = data_I_clean; data_Q = data_Q_clean;
M  = floor(length(data_I)/(Fs*Ts*Nf));

hasTrueRange=0; % only used with simulation

%% Simulation data
% outdir = sprintf('Results_simulations/Ts%.0fms_M%d_Nf%d_df%.0fMHz', Ts*1e3, M, Nf, df/1e6);
% datadir = sprintf('data_simulation/Ts%.0fms_M%d_Nf%d_df%.0fMHz', Ts*1e3, M, Nf, df/1e6);
% fname_re = sprintf('%s/cw_re_simulation.wav', datadir );
% fname_im = sprintf('%s/cw_im_simulation.wav', datadir );
% [data_I, Fs] = audioread(fname_re);
% [data_Q, ~]  = audioread(fname_im);

% % Try to load ground-truth data from the data generation file for comparison
% matdir = sprintf('%s/true_range.mat', datadir);
% hasTrueRange = exist(matdir, "file");
% if hasTrueRange
%     load(matdir);
% end

%% General code

R_amb = c / (2*df);
fprintf('Range ambiguity (resolution) = %.3f m\n', R_amb);
R_resol = c / (2*df*Nf);
fprintf('Range resolution = %.3f m\n', R_resol);

rx = complex(data_I, data_Q);
N  = length(rx);

samples_per_step  = round(Ts * Fs);
samples_per_sweep = Nf * samples_per_step;

% Reshape into sweep–frequency matrix
rx_matrix = reshape(rx(1:M*samples_per_sweep), samples_per_step, Nf, M);
S = squeeze(mean(rx_matrix, 1));  % [Nf x M]
freqs = f0 + (0:Nf-1)*df;

%% Range extraction (same as orignal processing but generalized for Nf > 2) via least-squares

% Compute phase of each tone
phases = angle(S);  % [Nf x M]

% Unwrap phases along frequency dimension
phases_unwrapped = unwrap(phases, [], 1);

% Estimate range via least-squares fit:
% phase = 4*pi*R*freq / c + const
R_est = zeros(1,M);
for m = 1:M
    p = phases_unwrapped(:,m);
    A = [4*pi*freqs(:)/c, ones(Nf,1)];
    x = A \ p;          % least squares: [R; phase_offset]
    R_est(m) = x(1);    % estimated range for this sweep (first value from the solution)
end

t_sweeps = ((0:M-1) * Ts);

%% Range–frequency FFT (range profile)
Nfft = length(S(1,:))*4;
range_axis = (0:Nfft-1) * (c / (2*Nfft*df));
range_profiles = fft(S, Nfft, 1);
range_magnitude = abs(range_profiles);
range_magnitude_norm = range_magnitude ./ max(range_magnitude);
range_magnitude_norm_dB = 20*log10(range_magnitude_norm);

figure;
imagesc(t_sweeps, range_axis, range_magnitude_norm_dB);
axis xy;
xlabel('Time (s)');
ylabel('Range (m)');
ylim([0 R_amb]);
title(sprintf('Range–Time Intensity (FFT over %d frequencies)', Nf));
colorbar;
colormap jet;
exportgraphics(gcf, sprintf('%s/RTI_FFT_Ts%.0fms_Nf%d_df%.0f.png', outdir, Ts*1e3, Nf, df/1e6), 'Resolution', 300);

%% Comparison with simulated true range / Estimated Range plot
figure;
if hasTrueRange
    subplot(2,1,1);
    plot(t_global, R_true, 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Range (m)');
    title('True Range'); grid on;

    subplot(2,1,2);
end
plot(t_sweeps, R_est, 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Estimated Range (m)');
title('Estimated Range'); grid on;

fname = sprintf('%s/%s_Ts%.0fms_Nf%d_df%.0fMHz.png', outdir, ternary(hasTrueRange, 'True_vs_Estimated_Range', 'Estimated_Range'), Ts*1e3, Nf, df/1e6);
exportgraphics(gcf, fname, 'Resolution', 300);


%% Extracts and unwraps the target range over time using peak tracking across sweeps

% get peak per sweep
[~, idx] = max(range_magnitude, [], 1);
R_mod = (idx-1) * (c / (2*Nfft*df));  % range in [0, Ramb)

% temporal unwrap using previous sweep to pick correct replica
Ramb = c/(2*df);
R_track = zeros(1,M);
R_track(1) = R_mod(1);
for m = 2:M
    k = round((R_track(m-1) - R_mod(m)) / Ramb);
    R_track(m) = R_mod(m) + k*Ramb;
end

figure;
if hasTrueRange
    subplot(2,1,1);
    plot(t_global, R_true, 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Range (m)');
    title('True Range'); grid on;
    subplot(2,1,2);
end

actual_performed_range = 2.15;
plot(t_sweeps, abs(R_track)./(max(abs(R_track))/actual_performed_range), 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Range (m)');
title('Unambiguous Tracked Range'); grid on;

fname = sprintf('%s/%s_Ts%.0fms_Nf%d_df%.0fMHz.png', outdir, ternary(hasTrueRange, 'True_vs_Tracked_Range', 'Tracked_Range'), Ts*1e3, Nf, df/1e6);
exportgraphics(gcf, fname, 'Resolution', 300);


%% Range–Time Intensity with R_track and R_est Tracking Line
% figure;
% imagesc(t_sweeps, range_axis, range_magnitude_norm_dB);
% hold on;
% plot(t_sweeps, R_track, 'c-', 'LineWidth',1.5);
% plot(t_sweeps, R_est, 'm-', 'LineWidth',1.5);
% axis xy;
% xlabel('Time (s)');
% ylabel('Range (m)');
% ylim([0 R_amb]);
% title(sprintf('Range–Time Intensity (FFT over %d frequencies)', Nf));
% colorbar;
% colormap jet;
% 
% exportgraphics(gcf, sprintf('%s/RTI_with_Rtrack_and_Rest_Ts%.0fms_Nf%d_df%.0f.png', outdir, Ts*1e3, Nf, df/1e6), 'Resolution', 300);

%% Range–Time Intensity with Max-Amplitude Tracking Line
% figure;
% imagesc(t_sweeps, range_axis, range_magnitude_norm_dB);
% axis xy;
% xlabel('Time (s)');
% ylabel('Range (m)');
% ylim([0 R_amb]);
% title(sprintf('Range–Time Intensity (FFT over %d frequencies)', Nf));
% colorbar;
% colormap jet;
% hold on;
% exportgraphics(gcf, sprintf('%s/RTI_with_MaxAmplitude_Ts%.0fms_Nf%d_df%.0f.png', outdir, Ts*1e3, Nf, df/1e6), 'Resolution', 300);
% 
% % Compute and overlay the maximum-amplitude line
% [~, idx_max] = max(range_magnitude_norm, [], 1);  % index of max for each time
% R_max = range_axis(idx_max);                      % convert indices to range
% plot(t_sweeps, R_max, 'w-', 'LineWidth', 1.5);    % white line
% 
% hold off;


%% Original processing: for Nf=2 steps
if Nf==2
    S1 = zeros(M,1);
    S2 = zeros(M,1);
    for m = 1:M
        base = (m-1)*samples_per_sweep + 1;
        idx1 = base : base + samples_per_step - 1;
        idx2 = base + samples_per_step : base + 2*samples_per_step - 1;
    
        if max(idx2) > length(rx)
            break
        end
    
        % Average IQ within each step to get one complex value per frequency
        S1(m) = mean(rx(idx1));
        S2(m) = mean(rx(idx2));
    
        % FFT-based dominant tone estimation
        % SIG1 = fft(rx(idx1));
        % SIG1=SIG1(1:length(SIG1)/2);
        % [~, idx_max] = max(abs(SIG1));        % find peak
        % S1(m) = SIG1(idx_max);              % store complex phasor
        % 
        % SIG2 = fft(rx(idx2));
        % SIG2=SIG2(1:length(SIG2)/2);
        % [~, idx_max] = max(abs(SIG2));
        % S2(m) = SIG2(idx_max);
        
    end
    
    dphi = angle(S2 .* conj(S1));          % (-pi, pi]
    
    R = dphi * (c / (4*pi*df));        % meters, wrapped phase
    R_unwrapped = unwrap(dphi) * (c / (4*pi*df));  % meters, unwrapped phase

    time_axis = (0:length(R)-1) * Ts;
    
    figure;
    plot(time_axis, abs(R_unwrapped)./(max(abs(R_unwrapped))/actual_performed_range), 'r-', 'LineWidth', 1);
    hold on;
    % plot(time_axis, R_unwrapped, 'r.-');
    % plot(time_axis, R, 'b.-');
    xlabel('Time (s)');
    ylabel('Unwrapped Range (m)');
    title('Unwrapped Range');
    grid on;
    hold off;
    exportgraphics(gcf, sprintf('%s/Unwrapped_Range_original_process_Ts%.0fms_Nf%d_df%.0f.png', outdir, Ts*1e3, Nf, df/1e6), 'Resolution', 300);

end


%% Helper function: used for inline conditional value selection (like an if-else in one line) to save the correct file name for both simulation and real cases
function out = ternary(cond, valTrue, valFalse)
    if cond
        out = valTrue;
    else
        out = valFalse;
    end
end


