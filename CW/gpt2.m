%% Clean
clear; clc;

%% Read audio (convert to .wav first if .m4a isn't supported by your MATLAB)
[data, Fs] = audioread("Range_Test_File2.wav");  % <-- use WAV if needed

%% Params
Tp         = 5e-3;                 % 5 ms up-chirp length (seconds)
f_start    = 2.405e9;
f_stop     = 2.495e9;
c          = 299792458;
N          = round(Tp * Fs);       % samples per up-chirp
delta_f    = f_stop - f_start;
delta_R    = c / (2*delta_f);
Rmax       = (N * delta_R) / 2;    % per slide
Fmax       = Fs/2;                 % not strictly used below

%% Channels
radar = data(:,1);       % backscatter
sync  = data(:,2);       % sync (square-like)

%% Robust sync thresholding + rising-edge detection
% Normalize sync and use hysteresis to avoid noise-triggered edges
s = sync - median(sync);                    % remove DC bias
hi = 0.2*max(abs(s)); lo = 0.1*max(abs(s));
is_hi = false(size(s));
state = false;
for k = 1:numel(s)
    if ~state && s(k) > hi, state = true;  end
    if  state && s(k) < lo, state = false; end
    is_hi(k) = state;
end

% Rising edges = start indices of up-chirps
rise_idx = find(diff([0; is_hi]) == 1);    % +1 transitions

%% Build MxN matrix: one row per complete up-chirp
rows = {};
for r = 1:numel(rise_idx)
    i0 = rise_idx(r);
    i1 = i0 + N - 1;
    if i1 <= numel(radar) && all(is_hi(i0:i1))  % stay within "up" region
        rows{end+1,1} = radar(i0:i1); %#ok<AGROW>
    end
end

if isempty(rows)
    error('No complete up-chirps found. Check Tp, Fs, or sync polarity.');
end

mat = cell2mat(rows.');   % size = M x N
M   = size(mat,1);

%% Mean Subtraction (column-wise per range cell), 2- and 3-pulse MTI
mat_ms   = mat - mean(mat, 1);                                   % Step 5
mat_2MTI = mat_ms(2:end,:) - mat_ms(1:end-1,:);
mat_3MTI = mat_ms(3:end,:) - 2*mat_ms(2:end-1,:) + mat_ms(1:end-2,:);

%% IFFT across columns with zero padding factor 4 (range transform)
ZP     = 4*N;
Y_raw  = ifft(mat,      ZP, 2);
Y_ms   = ifft(mat_ms,   ZP, 2);
Y_mti2 = ifft(mat_2MTI, ZP, 2);
Y_mti3 = ifft(mat_3MTI, ZP, 2);

% Convert to dB safely
to_dB  = @(X) 20*log10(abs(X)+eps);
Y_raw_dB  = to_dB(Y_raw);
Y_ms_dB   = to_dB(Y_ms);
Y_2_dB    = to_dB(Y_mti2);
Y_3_dB    = to_dB(Y_mti3);

%% Keep lower half in time (2N columns after 4N IFFT) per slide
Y_raw_dB  = Y_raw_dB(:,  1:2*N);
Y_ms_dB   = Y_ms_dB(:,   1:2*N);
Y_2_dB    = Y_2_dB(:,    1:2*N);
Y_3_dB    = Y_3_dB(:,    1:2*N);

%% Normalization (subtract matrix max) per slide
Y_raw_dB  = Y_raw_dB  - max(Y_raw_dB(:));
Y_ms_dB   = Y_ms_dB   - max(Y_ms_dB(:));
Y_2_dB    = Y_2_dB    - max(Y_2_dB(:));
Y_3_dB    = Y_3_dB    - max(Y_3_dB(:));

%% Axes
x_range = linspace(0, Rmax, 2*N);     % Step 11
t_axis  = (0:M-1)*Tp;                 % one time stamp per up-chirp row

%% Plots
figure(1); imagesc(x_range, t_axis, Y_raw_dB, [-50 0]);
xlabel('Range (m)'); ylabel('Time (s)'); title('Range: No MS, No MTI'); colorbar; xlim([0 100]);

figure(2); imagesc(x_range, t_axis, Y_ms_dB, [-50 0]);
xlabel('Range (m)'); ylabel('Time (s)'); title('Range: MS only'); colorbar; xlim([0 100]);

figure(3); imagesc(x_range, t_axis, Y_2_dB, [-50 0]);
xlabel('Range (m)'); ylabel('Time (s)'); title('Range: MS + 2-pulse MTI'); colorbar; xlim([0 100]);

figure(4); imagesc(x_range, t_axis, Y_3_dB, [-50 0]);
xlabel('Range (m)'); ylabel('Time (s)'); title('Range: MS + 3-pulse MTI'); colorbar; xlim([0 100]);
