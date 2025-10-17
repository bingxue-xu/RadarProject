%% SAR Time Domain Correlation Algorithm for Simulated Targets
% Group 7
% 23/09/2025
% reference: Lecture5 slide 43 - 46

%% Clearing
clear
clc

%% Radar values

fstart = 5.4e9;
fend = 5.5e9;
bw = fend - fstart;
Fs = 50e3;
Tp = 5e-3;
N = Tp*Fs;
fc = (fend + fstart)/2;
t = linspace(0, Tp, N);
c = 299792458;
lambda = c/fc;
spacing = lambda/2;
deltaR = c/(2*bw);
Rmax = deltaR*N/2;

%% Positions

x_pos = [-2050e-3 0e-3 2500e-3]; 
y_pos = [0 0 0]; 
z_pos = [10 50 90];

xp = -3000e-3:spacing:3000e-3;
yp = 0;

%% Acquiring simulated targets

sim_targ = SAR_Data_Generation(x_pos, y_pos, z_pos, xp, yp, N, t, fc, bw, Tp);

%% Creating range arrays

z_range = linspace(0, Rmax, N);
x_range = xp;
m = length(x_range);
n = length(z_range);

%% Evaluate reflectivity
f = zeros(m,n);
tic
for x_index = 1:length(x_range) %summation in x pixel
    for z_index = 1:length(z_range) %summation in z pixel
        x_img = x_range(x_index); % x pixel
        z_img = z_range(z_index); % z pixel
        pixel_value = 0;
        r = sqrt((x_img - x_range).^2 + z_img^2); % 1 x length(xp)
        t_d = 2 * r / c; % time delay of the received signal
        phase = exp(1i * 2 * pi * (fc * t_d.' + (bw / Tp) * t_d.' .* t));
        signal = squeeze(sim_targ(:, 1, :)); % 1xN matrix
        pixel_value = sum(sum(signal .* phase, 2)); % sum over t and x
        f(x_index, z_index) = abs(pixel_value);
    end
end

%% Multiplication by z_range^2

f = f.*(z_range.^2);
toc

%% Plotting

figure(1)
imagesc(z_range, x_range, abs(f/max(f(:))), [0 1]);
axis xy;
colorbar;
set(gca,'YDir','reverse');
xlabel('Z Down Range (m)');
ylabel('X Cross Range (m)');
title('SAR image, T_{p}=5ms, f_{start}=5.4GHz, f_{stop}=5.5GHz, F_{s}=50kHz');
saveas(gcf, 'images/SAR_fc.png');


