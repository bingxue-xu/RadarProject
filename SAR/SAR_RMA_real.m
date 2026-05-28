%% SAR Range Migration Algorithm (RMA) for Real Targets
% Group 7
% 26/09/2025
% reference: Lecture5 slide 54 - 68

%% Cleaning
clear
clc

% 1. load data
load("./radar_data/radar_4tables/radar_3tables.mat")
load("./radar_data/radar_4tables/sync_3tables.mat")
load("./radar_data/radar_4tables/Fs.mat")

clearvars -except radar_data sync_data Fs

%% Parameters
%% 3. set sync pulse to 1 if above 0
threshold = 0.01;
sync = sync_data > threshold;  % Convert sync to logical
Tp = 5e-3;                     % Pulse/chirp duration (s)
min_gap = round(Tp * Fs * 10); % Minimum number of consecutive zeros to count as a new position

N = round(Tp * Fs);    % expected number of samples per upchirp
if mod(N,2) ~= 0
    N = N+1;
end
tolerance = 20;              % tolerance

Nrp = round(0.5 * Fs);

f_start = 2.4e9;
f_end = 2.5e9;
BW = f_end - f_start;
f_c = (f_end + f_start) / 2;

spacing = 0.06;
c = 299792458;

%% 4. Robust grouping + protection against too-close starts
d = diff([0; sync(:); 0]);
starts_all = find(d == 1);
stops_all = find(d == -1) -1;
pulse_len_all = stops_all - starts_all +1;
valid_mask = abs(pulse_len_all - N) <= tolerance;

starts = starts_all(valid_mask);
stops = stops_all(valid_mask);

gaps = starts(2:end) - stops(1:end-1) -1;
gaps_threshold = 10*N;
separators = find(gaps > gaps_threshold);

% splite pulses into groups using separators
group_edges = [0, separators(:)', numel(gaps)];  % boundaries in the pulses index space
rail_positions = {};
rail_start_indices = [];

for k = 1:(length(group_edges)-1)
    first_pulse_idx = group_edges(k)+1;
    last_pulse_idx = group_edges(k+1)+1;
    pos_start = starts(first_pulse_idx);
    pos_stop = stops(last_pulse_idx);
    pos_len = pos_stop - pos_start +1;
    % require minimum position length
    if pos_len >= Nrp
        % clamp extremely long concatenation if you want: pos_stop = min(pos_start+max_pos_len-1, pos_stop);
        rail_positions{end+1} = radar_data(pos_start:pos_stop); %#ok<SAGROW>
        rail_start_indices(end+1) = pos_start; %#ok<SAGROW>
    end
end

%ensure odd number of positions as before
if mod(length(rail_positions), 2) == 0
    disp("Even number of rail positions detected. Dropping the last one.")
    rail_positions(end) = [];
    rail_start_indices(end) = [];
end

%% 5. Parsing values
radar_matrix = zeros(length(rail_start_indices), Nrp);
sync_matrix = zeros(length(rail_start_indices), Nrp);

%check for valid sync pulse
mean_thresh_min = 0.40;
mean_thresh_max = 0.60;
valid_rows = 1;

for i = 1:length(rail_start_indices)

    % % visualize sync and data alignment
    % idx = rail_start_indices(i)
    % window = 4*Nrp 
    % if (idx + window - 1) > length(radar_data)
    %     disp("skip: index out of range");
    %     continue
    % end
    % 
    % radar_segment = radar_data(idx:idx+window-1);
    % sync_segment = sync(idx:idx+window-1);
    % t = (0:length(radar_segment)-1)/Fs;
    % 
    % figure;
    % 
    % subplot(2,1,1);
    % plot(t, radar_segment, 'b');
    % title(sprintf("Radar backscatter, Rail %d (index %d)", i, idx));
    % xlabel("Time (s)");
    % ylabel("Amplitude");
    % grid on;
    % 
    % subplot(2,1,2);
    % plot(t, sync_segment, 'r');
    % title("Sync pulse");
    % xlabel("Time (s)");
    % ylabel("Logic value");
    % ylim([-0.1, 1.1]);
    % grid on;

    mean_sync = mean(sync(rail_start_indices(i)+Nrp:rail_start_indices(i) + 2*Nrp - 1));
    if mean_sync > mean_thresh_min && mean_sync < mean_thresh_max
        radar_matrix(valid_rows,:) = radar_data(rail_start_indices(i)+Nrp:rail_start_indices(i) + 2*Nrp - 1)';
        sync_matrix(valid_rows,:) = sync(rail_start_indices(i)+Nrp:rail_start_indices(i) + 2*Nrp - 1)';
        valid_rows = valid_rows + 1;
    end
end

radar_matrix = radar_matrix(1:valid_rows - 1, :);
sync_matrix = sync_matrix(1:valid_rows - 1, :);

%% RM4 3a. Plot the entire data and the sync signal (zoom in to illustrate) (Lec.5, slide 55).
N_samples = 1000*Fs;
t = (0:N_samples-1)/Fs;
% 
% figure;
% subplot(2,1,1);
% plot(t, radar_data(1:N_samples));
% title("Radar backscatter data");
% xlabel("Data Sample Number");
% ylabel("Amplitude");
% grid on;
% 
% subplot(2,1,2);
% plot(t, sync_data(1:N_samples));
% title("Sync data");
% xlabel("Data Sample Number");
% ylabel("Amplitude");
% grid on;
% saveas(gcf, 'images/sar_data.png');


%% try out on one rail position
i = 2;  
idx = rail_start_indices(i);
radar_segment = radar_data(idx : idx + 2*Nrp - 1);
sync_segment = sync(idx : idx + 2*Nrp - 1);

% find start idx of upchirp 
d_sync = diff([0; sync_segment; 0]); 
start_idx = find(d_sync == 1);
end_idx = find(d_sync == -1) - 1;

upchirps_this_pos = {};
for k = 1:numel(start_idx)
    if (start_idx(k) + N - 1) <= length(radar_segment)
        upchirps_this_pos{end+1} = radar_segment(start_idx(k):start_idx(k) + N - 1);
    end
end


%% processing upchirps 
% 1. parse all upchirps, 
% 2. average(integrate) them into one integrated-upchirp for each rail pos
% 3. stack integrated upchirps into a matrix -> data matrix 
%% Processing upchirps

upchirp_matrix = zeros(size(radar_matrix,1), N); %pre-allocate upchirp mat

% Fill upchirp_matrix
for i = 1:size(radar_matrix,1)
    radar_row = radar_matrix(i,:);
    sync_row  = sync_matrix(i,:);

    % Detect upchirp starts/stops
    d_sync = diff([0, sync_row, 0]);
    start_idx = find(d_sync == 1);
    end_idx   = find(d_sync == -1) - 1;

    upchirps_this_pos = zeros(numel(start_idx), N);
    counts = zeros(1,N);

    for k = 1:numel(start_idx)
        len_chirp = end_idx(k) - start_idx(k) + 1;
        if len_chirp >= N
            upchirps_this_pos(k,:) = radar_row(start_idx(k) : start_idx(k) + N - 1);
        else
            seg = radar_row(start_idx(k):end_idx(k));
            upchirps_this_pos(k,1:numel(seg)) = seg;
        end
    end

    upchirp_matrix(i,:) = mean(upchirps_this_pos, 1);

end

% figure;
% imagesc(abs(upchirp_matrix));
% xlabel('Sample index (fast time)');
% ylabel('Rail position (slow time)');
% title('Integrated Data Matrix (parsed by up-chirps)');
% colorbar;


%%
% samples_to_plot = min(1000, N);
% t = (0:samples_to_plot-1);
% 
% figure; hold on;
% for k = 1:length(upchirps_this_pos)
%     plot(upchirps_this_pos{k}(1:samples_to_plot));
% end
% title(sprintf("Multiple upchirps @ rail %d", i));
% xlabel("Sample index");
% ylabel("Amplitude");
% grid on;

% %% RM4 3b. 
% 
% Ns_to_plot = 2400;
% t = 0:Ns_to_plot-1;
% figure;
% 
% subplot(2,1,1);
% plot(t, sync_segment(1:Ns_to_plot));
% title("Position parsed sync data");
% xlabel("Data Sample Number (N_s)");
% ylabel("Amplitude");
% 
% subplot(2,1,2);
% plot(t, radar_segment(1:Ns_to_plot));
% title("Position parsed radar backscatter data");
% xlabel("Data Sample Number (N_s)");
% ylabel("Amplitude");
% 
% %% RM4 3c. 
% upchirp_matrix = cell2mat(upchirps_this_pos(:)).'; % rows = chirps
% 
% integrated_upchirp = mean(upchirp_matrix,1)
% figure;
% plot(integrated_upchirp);
% title(sprintf("Integrated up-chirp @ rail %d", i));
% xlabel("Sample index (N_s)");
% ylabel("Amplitude");
% grid on;


%% 6. Hilbert transform

fft_matrix = fft(upchirp_matrix,[],2);  %fft
data_matrix = ifft(fft_matrix, N/2, 2);  %ifft
data_matrix(isnan(data_matrix)) = 1e-30;

% %% RM4 3c. one rail postion before and after Hilbert transform
% % Pick one rail position, e.g. i = 10
% i = 10;
% 
% raw_chirp = upchirp_matrix(i,:);           % integrated up-chirp (real)
% % Hilbert transform → complex
% fft_chirp = fft(raw_chirp)
% fft_pos = fft_chirp(1:N/2);
% analytic_chirp = ifft(fft_pos, N/2);
% analytic_chirp(isnan(analytic_chirp)) = 1e-30;
% 
% t = 0:(N/2-1);  % sample index
% 
% figure;
% subplot(2,1,1);
% plot(0:N-1, raw_chirp);
% title(sprintf('Integrated up-chirp (before Hilbert) @ rail %d', i));
% xlabel('Data Sample Number (N_s)'); ylabel('Amplitude'); grid on;
% 
% subplot(2,1,2);
% plot(t, abs(analytic_chirp));
% title(sprintf('Integrated up-chirp (after Hilbert) @ rail %d', i));
% xlabel('Data Sample Number (N_s)'); ylabel('|Amplitude|'); grid on;



%% 7. Hann window

x = linspace(-N/4, N/4, N/2);
H = 0.5*(1 + cos(2*pi*x*(2/N)));
data_hann = data_matrix .* H; 

% %% RM4 3e. Show results with and without using the Hann window.
i = 10;
t = 0:(N/2-1);

% figure;
% subplot(1,2,1);
% plot(t, data_matrix(i,:));
% title('Data from one row (before hann windowing)');
% xlabel('Data Sample Number (N_s)');
% ylabel('Amplitude');
% grid on;
% 
% subplot(1,2,2);
% plot(t, data_hann(i,:));
% title('Data from one row (after hann windowoing)');
% xlabel('Data Sample Number (N_s)');
% ylabel('Amplitude');
% grid on;
% saveas(gcf, 'images/SAR_hann_plot.png');
% 
% 
% % Compare full data matrices (imagesc)
% figure;
% subplot(1,2,1);
% imagesc(abs(data_matrix));
% title('Data Matrix (no Hann window)');
% xlabel('Fast-time samples (N_s)');
% ylabel('Rail position'); colorbar;
% 
% subplot(1,2,2);
% imagesc(abs(data_hann));
% title('Data Matrix (with Hann window)');
% xlabel('Fast-time samples (N_s)');
% ylabel('Rail position'); colorbar;
% saveas(gcf, 'images/SAR_hann_matrices.png');


%% 8. Setting parameters
t = linspace(0, Tp, N/2);
L = spacing * size(data_hann, 1);
x_p = linspace(-L/2, L/2, L/spacing);

%% 9. Generate kt
k_t = (2*pi*f_c)/c + (2*pi*BW*t)/(Tp*c);

% %% 10. plot check
% figure()
% imagesc(k_t, x_p, angle(data_matrix), [-3 3]);
% axis xy;
% colorbar;
% set(gca,'YDir','reverse');
% xlabel('k_t (rad/m)');
% ylabel('Synthetic Aperature Position, x_p (m)');
% title('Phase Before Along Track FFT');

%% 11. zero padding
num_zero = 2049;
zeros_matrix = zeros(num_zero, N/2);
data2_matrix = [zeros_matrix;data_hann;zeros_matrix];
%% 12. update L and xp
L = spacing * size(data2_matrix, 1);
x_p = linspace(-L/2, L/2, L/spacing);

%% 13-16. k_x and 2-D matrix, k_t 2-D matrix, kz
delta_k_x = (2*pi)/(2*max(x_p));
k_x = delta_k_x * (-(size(x_p,2) - 1) / 2 : (size(x_p,2) - 1 )/ 2);

k_x_2d = repmat(k_x.', [1,N/2]);
k_t_2d = repmat(k_t, [size(x_p,2),1]);
k_z_2d = (4 * (k_t_2d .^ 2) - (k_x_2d .^ 2)) .^ 0.5;

%% 17. Stolt interpolation
k_z_1d_uni = linspace(min(min(k_z_2d)), max(max(k_z_2d)), size(k_z_2d, 2));
k_z_2d_uni = repmat (k_z_1d_uni, [size(x_p,2),1]);

% %% 18. Plot k-space
% figure()
% scatter(k_x_2d(1:10:end), k_z_2d(1:10:end), 1, 'Marker','.'); hold on
% scatter(k_x_2d(1:10:end), k_z_2d_uni(1:10:end), 1, 'Marker','.'); hold off

%% 19. FFT
%plot magnitude in dB and phase
S_B = fftshift(fft(data2_matrix,[],1),1);
% %% 20. Plot magnitude and phase
% figure()
% imagesc(k_t, k_x, 20*log10(abs(S_B)), [-15 20]);
% colorbar;
% figure()
% imagesc(k_t, k_x, angle(S_B), [-3 3]);
% colorbar;

%% 21. Stolt interpolation
S_B2 = zeros(size(S_B,1),N/2);
for i = 1:size(S_B, 1)
    S_B2(i,:) = interp1(k_z_2d(i,:), S_B(i,:), k_z_1d_uni);
end
S_B2(isnan(S_B2)) = 1e-30;

% %% 22. Plot phase after stolt
% figure()
% imagesc(k_z_1d_uni, k_x, angle(S_B2));
% colorbar;

%% 23. 2-D ifft
ft_1 = ifft(S_B2,[],1);
ft_2 = ifft(ft_1,[],2);

%% 24. Flip and rotate
data3_matrix = fliplr(rot90(ft_2, 1));

%% Indexes to trunc
delta_fz = c*(max(k_z_1d_uni) - min(k_z_1d_uni))/(4*pi);
delta_k_z = max(k_z_1d_uni) - min(k_z_1d_uni);

Rmax = N*pi/delta_k_z;
Rail_rmax = size(ft_2,1)*spacing;

d_range_1 = 1;
d_range_2 = floor(Rmax/2);
c_range_1 = -10;
c_range_2 = 10;

d_index_1 = ceil((size(data3_matrix,1)/Rmax)*d_range_1);
d_index_2 = ceil((size(data3_matrix,1)/Rmax)*d_range_2);

c_index_1 = ceil((size(data3_matrix,2)/Rail_rmax)*(c_range_1 + (Rail_rmax/2)));
c_index_2 = ceil((size(data3_matrix,2)/Rail_rmax)*(c_range_2 + (Rail_rmax/2)));

truncated_data_matrix = data3_matrix(d_index_1:d_index_2,c_index_1:c_index_2);

%% 25. Downrange and crossrange vector
downrange = linspace(-1*d_range_1, -1*d_range_2, size(truncated_data_matrix, 1));
crossrange = linspace(c_range_1, c_range_2, size(truncated_data_matrix, 2));

%% 26. Column multiplication
truncated2_data_matrix = truncated_data_matrix .* ((downrange.^2)');

%% 27. dB conversion and plot
image = 20*log10(abs(truncated2_data_matrix));
figure()
imagesc(crossrange, downrange, image, [max(max(image))-25, max(max(image))+0]);
colorbar;
set(gca,'YDir','reverse');
xlabel('Crossrange (meter)');
ylabel('Downrange (meter)');
title('SAR real measurement final image');
saveas(gcf, 'images/sar_final_image.png');

%% RM4 3f. Show results with and without multiplication of the down range data with the square of the down range vector.
% image_no_comp = 20*log10(abs(truncated_data_matrix));
% image = 20*log10(abs(truncated2_data_matrix));
% figure()
% subplot(1,2,1);
% imagesc(crossrange, downrange, image_no_comp, [max(max(image_no_comp))-25, max(max(image_no_comp))+0]);
% colorbar;
% set(gca,'YDir','reverse');
% xlabel('Crossrange (meter)');
% ylabel('Downrange (meter)');
% title('Final Image no R^2 compensation');
% saveas(gcf, 'images/SAR_REAL_no_R².png');
% 
% subplot(1,2,2);
% imagesc(crossrange, downrange, image, [max(max(image))-25, max(max(image))+0]);
% colorbar;
% set(gca,'YDir','reverse');
% xlabel('Crossrange (meter)');
% ylabel('Downrange (meter)');
% title('Final Image with R^2 compensation');
% saveas(gcf, 'images/SAR_REAL_R².png');
