%Function to create targets for SAR processing

%Input position of the targets in the azimuth, elevation and range space.
%x_pos, y_pos and z_pos as vectors e.g. for two targets at 80cm and 30cm 
%in x and 0cm and 0cm in y and 1m and 5m in range write x_pos=[80e-2 30e-2], 
% y_pos=[0 0], z_pos=[1 5]

%Input radar positions on a rail in azimuth and elevation.
%xp is vector for radar x positions with spacing of half wavelength between
%each postion so movement of 2m would be written as xp=-1:lambda/2:1. 
%yp is the vector for radar y position. Set the yp=0 as we are not moving 
%in elevation.

%N is number of samples N = Tp*Fs

%t is time array t = linspace(0, Tp, N)

%fc is the center frequency

%BW is teh bandwidth

%Tp is the pulse width or the up-chirp time

function SAR_Data_Time_Domain = SAR_Data_Generation(x_pos, y_pos, z_pos, xp, yp, N, t, fc, BW, Tp)

c = 299792458; %(m/s) speed of light
f = 1*ones(1, length(x_pos)); % reflectivity function for the targets

n_targets = length(x_pos); % Total number of targets
Tot_x = length(xp); % Total number of x positions/measurements
Tot_y = length(yp); % Total number of y positions/measurements

SAR_Data_Time_Domain = zeros(Tot_x,Tot_y,N);

for nt = 1:n_targets
    for m = 1:length(xp)
        r = sqrt((x_pos(nt)-xp(m)).^2 + y_pos(nt).^2 + z_pos(nt).^2);
        tau = 2 * r / c;
        phase = cos((2*2*pi*fc*r/c)+(2*2*pi*BW*r*t/(c*Tp)));  % 1×N
        SAR_Data_Time_Domain(m, 1, :) = SAR_Data_Time_Domain(m, 1, :) + reshape(f(nt)*(1./(r.^2)).*phase, 1, 1, []);
    end
end