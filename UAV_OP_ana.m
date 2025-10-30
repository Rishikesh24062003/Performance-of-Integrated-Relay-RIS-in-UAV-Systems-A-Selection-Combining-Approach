%% Clear and close all
clear;
close all;
clc;

%% Parameters
M = 1e6; % Number of samples (for simulation)
Average_SNR_dB = -30:5:50; % Average SNR in dB
Average_SNR = 10.^(Average_SNR_dB ./ 10); % Convert SNR to linear scale
Threshold_SNR_dB = -30; % Threshold SNR in dB
Gamma_0 = 10^(Threshold_SNR_dB / 10); % Convert threshold SNR to linear scale

% Nakagami-m parameters
m = 2; % Shape parameter
Omega_UAV = 1; % Mean SNR of UAV link
Omega_IRS = 1; % Mean SNR of IRS link

% Gains
G_s = 10; G_u = 10; G_r = 10; G_d = 10;

% Path Losses
d1 = 10; % Distance from source to UAV in meters
d2 = 10; % Distance from UAV to destination in meters
f = 1.6e9; % Frequency in Hz (1.6 GHz)
c = 3e8; % Speed of light in m/s
lambda = c / f; % Wavelength

% Free Space Path Loss parameters
PL_UAV = (G_s * G_u * lambda^2) / (4 * pi * d1)^2; % Source to UAV path loss
PL_IRS = (G_r * G_d * lambda^2) / (4 * pi * d2)^2; % UAV to Destination path loss

% Initialize arrays for analytical outage probability results
OutageProb_UAV = zeros(size(Average_SNR_dB));
OutageProb_IRS = zeros(size(Average_SNR_dB));
OutageProb_combined = zeros(size(Average_SNR_dB));

%% Analytical outage probability calculation loop for UAV and IRS

for jj = 1:length(Average_SNR_dB)
    % Analytical Outage Probability Calculation for UAV
    gamma_th_uav = Gamma_0 / (Average_SNR(jj) * PL_UAV); % Normalized threshold for Nakagami-m CDF
    F_gamma_UAV = gammainc(m * gamma_th_uav, m, 'lower');
    OutageProb_UAV(jj) = F_gamma_UAV;
    
    % Analytical Outage Probability Calculation for IRS
    gamma_th_irs = Gamma_0 / (Average_SNR(jj) * PL_IRS); % Normalized threshold for Nakagami-m CDF
    F_gamma_IRS = gammainc(m * gamma_th_irs, m, 'lower');
    OutageProb_IRS(jj) = F_gamma_IRS;
    
    % Combined Outage Probability using Selection Combining
    OutageProb_combined(jj) = OutageProb_UAV(jj) * OutageProb_IRS(jj);
end

%% Plotting Outage Probability results
figure;
semilogy(Average_SNR_dB, OutageProb_combined, 'k--', 'LineWidth', 1.5);
xlabel('Average SNR (dB)');
ylabel('Outage Probability');
title('Outage Probability vs Average SNR for UAV Relay and IRS with Selection Combining');
legend('Combined');
grid on;
axis([-30 50 1e-6 1]);
