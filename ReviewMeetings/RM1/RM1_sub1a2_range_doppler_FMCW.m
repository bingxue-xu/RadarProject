%%clean up
clear; 
clc;

%%params
c = 3e8; % speed of light (m/s)
f_start = 2.405e9; % start freq (Hz)
f_stop = 2.495e9; % stop freq (Hz)
deltaf = f_stop - f_start; % sweep bandwidth (Hz)
Tp = 5e-3; % chirp duration (s) 
fc = 0.5*(f_start + f_stop); % center frequency (Hz)
lambda = c / fc; % wavelength (m)

%%read and sep data
filename = "Range_Test_File2.m4a";
[data, Fs] = audioread(filename);
radar_backscatter = data(:,1);
sync_signal = data(:,2);

chunk = Fs * 1;
Indx = 1;
threshold = 0;

%%video output stuff
v = VideoWriter('RangeDopplerSpectrogram.avi');
v.FrameRate = 5;
open(v);
fig = figure('Color','w');

%%loop
while Indx + chunk - 1 <= length(radar_backscatter)

    % 1-sec chunks
    rb_chunk = radar_backscatter(Indx : Indx + chunk - 1);
    sync_chunk = sync_signal(Indx : Indx + chunk - 1);

    syncPulse = zeros(size(sync_chunk));
    syncPulse(sync_chunk > threshold) = 1;

    d_sync = diff([0; syncPulse; 0]); %add zeros incase mid up/down chirp for detection
    start_idx = find(d_sync == 1);
    end_idx = find(d_sync == -1) - 1;

    % Make all chirps same len
    N_chirp = round(Fs * Tp); % fixed chirp length
    num_chirps = numel(start_idx);

    upchirp_matrix = zeros(num_chirps, N_chirp);
    for k = 1:num_chirps
        %zero pad short segnments
        if end_idx(k) - start_idx(k) + 1 >= N_chirp
            upchirp_matrix(k,:) = rb_chunk(start_idx(k) : start_idx(k) + N_chirp - 1);
        else
            %if a chirp is shorter than expected -> 0-pad
            seg = rb_chunk(start_idx(k):end_idx(k));
            upchirp_matrix(k,1:numel(seg)) = seg;
        end
    end

    %MS clutter rejection
    upchirp_ms = upchirp_matrix - mean(upchirp_matrix, 1);

    %2 pulse MTI
    data_2MTI = upchirp_ms(2:end, :) - upchirp_ms(1:end-1, :);

    %IFFT's
    Y1 = ifft(data_2MTI, [], 2); %fast-time IFFT so range
    Y2 = ifft(Y1, [], 1); %slow-time IFFT for vel
    Y2abs = abs(Y2);
    Y2_dB = 20*log10(Y2abs + eps);

    %kepp lower half and normalize step
    N = size(Y2_dB, 2);          
    Y2_lower = Y2_dB(:, 1:floor(N/2));          
    Y2_lower = Y2_lower - max(Y2_lower(:));          

    %vel axis stuff
    d_f = linspace(0, 1/(2*Tp), size(Y2_lower,1)); %Doppler frequency (Hz)
    velocity = (lambda/2) * d_f;

    %range using deltaR and Rmax
    deltaR = c / (2 * deltaf);
    Rmax = (N * deltaR) / 2;
    range = linspace(0, Rmax, floor(N/2));          

    % plot frame
    imagesc(range, velocity, Y2_lower, [-40 0]);  
    axis ij;
    xlabel('Range (m)'); ylabel('Velocity (m/s)');
    title('Range–Doppler Plot (IFFT–IFFT, MS + 2-pulse MTI)');
    colorbar;
    xlim([0 100]);

    drawnow;
    frame = getframe(fig);
    writeVideo(v, frame);

    %indexing obv
    Indx = Indx + chunk;
end

close(v);
disp('Saved video: RangeDopplerSpectrogram.avi');
