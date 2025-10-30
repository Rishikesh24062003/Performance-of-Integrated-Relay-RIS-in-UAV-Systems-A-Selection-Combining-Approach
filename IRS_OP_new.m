%% Clear and close all
% clear;
close all;
clc;

%% Parameters for Nakagami-m fading channel with N RIS elements
N = 16; % Number of reflecting elements
M = 1e6; % Number of samples
Average_SNR_dB = -10:5:50; % Average SNR in dB
Average_SNR = 10.^(Average_SNR_dB./10); % Convert SNR to linear scale
Threshold_SNR_dB = -10; % Threshold SNR in dB
Threshold_SNR = 10^(Threshold_SNR_dB/10); % Convert threshold SNR to linear scale

m = 2; % Shape parameter for Nakagami-m distribution
omega = 2; % Scale parameter for Nakagami-m distribution

% Carrier frequency and path-loss calculation
fc = 3e9; % Carrier frequency (3 GHz)
c = 3e8; % Speed of light
lambda = c / fc; % Wavelength
p = lambda / (4 * pi); % Path loss factor

% Antenna gains
antenna_gain_S = db2pow(10); % Source antenna gain in linear scale
antenna_gain_RIS = db2pow(10); % RIS element gain in linear scale
antenna_gain_D = db2pow(10); % Destination antenna gain in linear scale

d1 = 10; % Distance from source to RIS
d2 = 10; % Distance from RIS to destination

% Path-loss calculation
path_loss_h = (antenna_gain_S * antenna_gain_RIS * p^2 / d1^2)^-1;
path_loss_g = (antenna_gain_RIS * antenna_gain_D * p^2 / d2^2)^-1;
path_loss = path_loss_h * path_loss_g;

% Analytical outage probability calculation
theta = sqrt((m * m) / (omega * omega));
E_vk1 = (gamma(m + 0.5) * gamma(m + 0.5)) * theta^(-1)/ (gamma(m) * gamma(m));
E_vk2 = (gamma(m + 1) * gamma(m + 1)) * theta^(-2)/ (gamma(m) * gamma(m));
var_vk = E_vk2 - E_vk1^2;
E_zk = N * E_vk1;
var_zk = N * var_vk;
ak = (E_zk^2) / var_zk;
bk = var_zk / E_zk;

Pout_RF = zeros(size(Average_SNR));
for ii = 1:length(Average_SNR)
    sk = sqrt((path_loss * Threshold_SNR) / (Average_SNR(ii) * bk^2));
    Pout_RF(ii) = gammainc(sk, ak, 'lower'); % Analytical outage probability
end

% Simulation of outage probability
hi = sqrt(gamrnd(m, omega/m, N, M)); % Nakagami-m fading channel coefficients (hi)
gi = sqrt(gamrnd(m, omega/m, N, M)); % Nakagami-m fading channel coefficients (gi)
h1 = (1 / sqrt(path_loss)) * hi .* gi; % Nakagami-m fading channel (h1)

h = sum(h1); % Sum of channel effects over N elements

P_out = zeros(size(Average_SNR_dB));
for k = 1:length(Average_SNR_dB)
    instt_SNR = Average_SNR(k) * (abs(h).^2); % Instantaneous SNR
    count = sum(instt_SNR < Threshold_SNR); % Count of SNR below threshold
    P_out(k) = count / M; % Simulated outage probability
end

% Plotting the results
figure;
semilogy(Average_SNR_dB, P_out, 'g-o', 'LineWidth', 1.5); % Simulation
hold on;
grid on;
semilogy(Average_SNR_dB, Pout_RF, 'k--', 'LineWidth', 1.5); % Analytical
xlabel('Average SNR (dB)');
ylabel('Outage Probability');
title('Outage Probability vs Average SNR for Nakagami-m Fading Channel with IRS');
legend( 'Simulation','Analytical');
grid on;
hold off;
axis([10 60 1e-6 1])