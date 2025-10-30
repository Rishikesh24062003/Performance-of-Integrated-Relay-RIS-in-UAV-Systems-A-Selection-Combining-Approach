%% Clear and close all
clear;
close all;
clc;

%% Parameters for Relay System
M = 1e6; % Number of samples (for simulation)
Average_SNR_dB = -30:5:50; % Average SNR in dB
Average_SNR = 10.^(Average_SNR_dB ./ 10); % Convert SNR to linear scale
Threshold_SNR_dB = -30; % Threshold SNR in dB
Gamma_0 = 10^(Threshold_SNR_dB / 10); % Convert threshold SNR to linear scale

% Nakagami-m parameters
m = 2; % Shape parameter
Omega_SR = 1; % Mean SNR of Source-Relay link
Omega_RD = 1; % Mean SNR of Relay-Destination link

% Gains
G_s = 10; G_r = 10; G_d = 10;

% Path Losses
d1 = 10; % Distance from source to relay in meters
d2 = 10; % Distance from relay to destination in meters
f = 1.6e9; % Frequency in Hz (1.6 GHz)
c = 3e8; % Speed of light in m/s
lambda = c / f; % Wavelength

% Free Space Path Loss parameters
PL1 = (G_s * G_r * lambda^2) / (4 * pi * d1)^2; % Source to Relay path loss
PL2 = (G_r * G_d * lambda^2) / (4 * pi * d2)^2; % Relay to Destination path loss

% Initialize arrays for outage probability simulation and analytical results
OutageProb_sim = zeros(size(Average_SNR_dB ));
OutageProb_analytical = zeros(size(Average_SNR_dB ));

%% Simulation loop over different transmit powers
for jj = 1:length(Average_SNR_dB )
    % Generate Nakagami-m fading coefficients
    h1 = sqrt(gamrnd(m, 1/m, 1, M)); % Nakagami-m fading coefficients for h1
    h2 = sqrt(gamrnd(m, 1/m, 1, M)); % Nakagami-m fading coefficients for h2

    % Apply path losses including FSPL
    h1f = h1;
    h2f = h2;

    % Effective SNR considering path losses and transmit power
    SNR1 = Average_SNR(jj) * h1f.^2 *(sqrt(PL1));
    SNR2 = Average_SNR(jj) * h2f.^2*(sqrt(PL2))  ;
    SNR = min(SNR1, SNR2); % Minimum SNR for outage calculation

    % Calculate outage probability from simulation
    OutageProb_sim(jj) = mean(SNR < Gamma_0);
end

%% Analytical outage probability calculation loop

for jj = 1:length(Average_SNR_dB )
    % Compute CDF for Source-Relay link
% Analytical Outage Probability Calculation
    gamma_th_sr =  Gamma_0  /Omega_SR;  % Normalized threshold for Nakagami-m CDF
    F_gamma_SR = 1-gammainc(m *(1/sqrt(PL1))* gamma_th_sr./Average_SNR, m, 'lower');

    gamma_th_rd =  Gamma_0  /Omega_RD;  % Normalized threshold for Nakagami-m CDF
    F_gamma_RD = 1-gammainc(m *(1/sqrt(PL2))* gamma_th_rd./Average_SNR, m, 'lower');

    % Compute outage probability analytically
    OutageProb_analytical = 1-(F_gamma_SR.*F_gamma_RD);
end

%% Plotting Outage Probability results
figure;
semilogy(Average_SNR_dB , OutageProb_sim, 'r-s', 'LineWidth', 1.5);
hold on;
semilogy(Average_SNR_dB , OutageProb_analytical, 'k--', 'LineWidth', 1.5);
xlabel('Transmit Power (dB)');
ylabel('Outage Probability');
title('Outage Probability vs Transmit Power for Nakagami-m Fading Channel with Relay System (Including FSPL)');
legend('Simulation', 'Analytical');
grid on;
axis([-30 10 1e-6 1])