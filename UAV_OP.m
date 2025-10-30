%% Clear and close all
clear;
close all;
clc;

%% Parameters
N = 16; % Number of reflecting elements
M = 1e6; % Number of samples
Average_SNR_dB = -30:5:50; % Average SNR in dB
Average_SNR = 10.^(Average_SNR_dB ./ 10); % Convert SNR to linear scale
Threshold_SNR_dB = -10; % Threshold SNR in dB
Threshold_SNR = 10^(Threshold_SNR_dB / 10); % Convert threshold SNR to linear scale

m = 2; % Shape parameter for Nakagami-m distribution
omega = 2; % Scale parameter for Nakagami-m distribution

% Carrier frequency and path-loss calculation
fc = 3e9; % Carrier frequency (3 GHz)
c = 3e8; % Speed of light
lambda = c / fc; % Wavelength
p = lambda / (4 * pi); % Path loss factor

%% Parameters for IRS System
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

%% Analytical Outage Probability Calculation for IRS
theta = sqrt((m * m) / (omega * omega));
E_vk1 = (gamma(m + 0.5) * gamma(m + 0.5)) * theta^(-1) / (gamma(m) * gamma(m));
E_vk2 = (gamma(m + 1) * gamma(m + 1)) * theta^(-2) / (gamma(m) * gamma(m));
var_vk = E_vk2 - E_vk1^2;
E_zk = N * E_vk1;
var_zk = N * var_vk;
ak = (E_zk^2) / var_zk;
bk = var_zk / E_zk;

Pout_IRS = zeros(size(Average_SNR));
for ii = 1:length(Average_SNR)
    sk = sqrt((path_loss * Threshold_SNR) / (Average_SNR(ii) * bk^2));
    Pout_IRS(ii) = gammainc(sk, ak, 'lower'); % Analytical outage probability for IRS
end

%% Simulation of Outage Probability for IRS
hi = sqrt(gamrnd(m, omega/m, N, M)); % Nakagami-m fading channel coefficients (hi)
gi = sqrt(gamrnd(m, omega/m, N, M)); % Nakagami-m fading channel coefficients (gi)
h1 =  hi .* gi; % Nakagami-m fading channel (h1)

h = sum(h1); % Sum of channel effects over N elements

P_out_IRS_sim = zeros(size(Average_SNR_dB));
for k = 1:length(Average_SNR_dB)
    instt_SNR = Average_SNR(k) * (abs(h).^2); % Instantaneous SNR
    count = sum(instt_SNR < Threshold_SNR); % Count of SNR below threshold
    P_out_IRS_sim(k) = count / M; % Simulated outage probability for IRS
end

%% Parameters for Relay System
Average_SNR_dB_Relay = Average_SNR_dB; % Average SNR in dB for relay system
Average_SNR_Relay = 10.^(Average_SNR_dB_Relay ./ 10); % Convert SNR to linear scale
Threshold_SNR_Relay = 10^(Threshold_SNR_dB / 10); % Convert threshold SNR to linear scale

% Path Losses
PL1 = (antenna_gain_S * antenna_gain_RIS * lambda^2) / (4 * pi * d1)^2; % Source to Relay path loss
PL2 = (antenna_gain_RIS * antenna_gain_D * lambda^2) / (4 * pi * d2)^2; % Relay to Destination path loss

%% Simulation loop over different transmit powers for Relay System
OutageProb_sim_Relay = zeros(size(Average_SNR_dB_Relay));
for jj = 1:length(Average_SNR_dB_Relay)
    % Generate Nakagami-m fading coefficients for Relay system
    h1 = sqrt(gamrnd(m, omega/m, 1, M)); % Nakagami-m fading coefficients for h1 (Source-Relay link)
    h2 = sqrt(gamrnd(m, omega/m, 1, M)); % Nakagami-m fading coefficients for h2 (Relay-Destination link)

    % Apply path losses including Free-Space Path Loss (FSPL)
    h1f = h1 * PL1;
    h2f = h2 * PL2;

    % Effective SNR considering path losses and transmit power
    SNR1 = Average_SNR_Relay(jj) * h1f.^2; % SNR for Source-Relay link
    SNR2 = Average_SNR_Relay(jj) * h2f.^2; % SNR for Relay-Destination link
    SNR_Relay = min(SNR1, SNR2); % Minimum SNR for outage calculation

    % Calculate outage probability from simulation for Relay system
    OutageProb_sim_Relay(jj) = mean(SNR_Relay < Threshold_SNR_Relay);
end

%% Analytical Outage Probability Calculation for Relay System
OutageProb_analytical_Relay = zeros(size(Average_SNR_dB_Relay));
for jj = 1:length(Average_SNR_dB_Relay)
    % Compute CDF for Source-Relay link
    gamma_th_sr = Threshold_SNR_Relay / omega;
    F_gamma_SR = 1 - gammainc(m * gamma_th_sr / Average_SNR_Relay(jj), m, 'lower');

    % Compute CDF for Relay-Destination link
    gamma_th_rd = Threshold_SNR_Relay / omega;
    F_gamma_RD = 1 - gammainc(m * gamma_th_rd / Average_SNR_Relay(jj), m, 'lower');

    % Compute outage probability analytically for Relay system
    OutageProb_analytical_Relay(jj) = 1 - (F_gamma_SR * F_gamma_RD);
end

%% Selection Combining: Minimum of IRS and Relay outage probabilities
Pout_SC_sim = P_out_IRS_sim .* OutageProb_sim_Relay; % Selection Combining Simulation

Pout_SC_analytical = Pout_IRS .* OutageProb_analytical_Relay; % Selection Combining Analytical

%% Plotting Selection Combining Outage Probability results
figure;
semilogy(Average_SNR_dB, Pout_SC_sim, 'm-^', 'LineWidth', 1.5); % Selection Combining Simulation
hold on;
semilogy(Average_SNR_dB, Pout_SC_analytical, 'c--', 'LineWidth', 1.5); % Selection Combining Analytical
xlabel('Average SNR (dB)');
ylabel('Outage Probability');
title('Outage Probability vs Average SNR for Integrated UAV-IRS Relaying with Selection Combining');
legend('Selection Combining Simulation', 'Selection Combining Analytical');
grid on;
axis([-30 50 1e-6 1]);
hold off;

% Display the plot
show();
