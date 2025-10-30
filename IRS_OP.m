clear all;
clc;

% Parameters
M = 100000; % Reduced number of symbols for testing
Pt_dB = -30:5:50; % Transmit power in dB
Pt = 10.^(Pt_dB/10); % Transmit power in linear scale
No = 1; % Noise power
m = 2; % Nakagami-m fading parameter
omega = 1; % Omega parameter for Nakagami-m distribution

% Threshold for outage probability
SNR_th_dB = -100; % Threshold SNR in dB
SNR_th = 10^(SNR_th_dB/10);

% Free Space Path Loss parameters
d1 = 10; % Distance from source to IRS in meters
d2 = 20; % Distance from IRS to destination in meters
f = 1.6e9; % Frequency in Hz (6 GHz)
c = 3e8; % Speed of light in m/s
G_S = 1; % Gain of source antenna
G_D = 1; % Gain of destination antenna

% Calculate wavelength
lambda = c / f;

% Calculate Path Loss for IRS
PL_IRS = (lambda^2 * sqrt(G_S * G_D)) / (16 * pi * d1 * d2);

% Generate BPSK symbols
ip = rand(1, M) > 0.5; % Generate 0,1 with equal probability
s = 2 * ip - 1; % BPSK modulation 0 -> -1, 1 -> 1 

% Initialize the outage probability array
OutageProb = zeros(1, length(Pt_dB));

% IRS Simulation
N = 32; % Number of reflecting elements

% Generate Nakagami-m fading coefficients for N elements (outside loop)
hi = sqrt(gamrnd(m, omega/m, [N, M]));  % |h|^2 is gamma distributed, |h| is sqrt(gamma)
gi = sqrt(gamrnd(m, omega/m, [N, M]));  % |g|^2 is gamma distributed, |g| is sqrt(gamma)

% Calculate the combined channel effect (outside loop)
h = PL_IRS * sum(hi .* gi);

for jj = 1:length(Pt)
    % Calculate instantaneous SNR
    instt_SNR = Pt(jj) * (abs(h).^2) / No;

    % Calculate outage probability
    OutageProb(jj) = mean(instt_SNR < SNR_th);
end

% Plotting Outage Probability results
figure;
semilogy(Pt_dB, OutageProb, 'r-s', 'LineWidth', 1.5);
xlabel('Transmit Power (dB)');
ylabel('Outage Probability');
title('Outage Probability vs Transmit Power for Nakagami-m Fading Channel with IRS');
grid on;
