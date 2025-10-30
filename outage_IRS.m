clear all;
clc;

% Parameters
M = 10^6; % Number of symbols
Pt_dB = -20:2:30; % Transmit power in dB
Pt = 10.^(Pt_dB/10); % Transmit power in linear scale
No = 1; % Noise power
m = 2; % Nakagami-m fading parameter
omega = 1; % Omega parameter for Nakagami-m distribution
N = 32; % Number of reflecting elements in IRS
SNR_th_dB = 10; % Threshold SNR in dB
SNR_th = 10^(SNR_th_dB/10); % Threshold SNR in linear scale

% Initialize the outage probability array
OutageProb = zeros(1, length(Pt_dB));

for jj = 1:length(Pt) 
    % Generate Nakagami-m fading coefficients for IRS with N elements
    hi = sqrt(gamrnd(m, omega/m, [N, M]));  % |h|^2 is gamma distributed, |h| is sqrt(gamma)
    gi = sqrt(gamrnd(m, omega/m, [N, M]));  % |g|^2 is gamma distributed, |g| is sqrt(gamma)
    h_IRS = sum(hi .* gi); % Combined IRS channel

    % Calculate instantaneous SNR
    inst_SNR = Pt(jj) * (abs(h_IRS).^2) / No;

    % Calculate outage probability
    OutageProb(jj) = mean(inst_SNR < SNR_th);
end

% Plotting Outage Probability results
figure;
semilogy(Pt_dB, OutageProb, 'r-s', 'LineWidth', 1.5);
xlabel('Transmit Power (dB)');
ylabel('Outage Probability');
title('Outage Probability vs Transmit Power for Nakagami-m Fading Channel with IRS');
grid on;