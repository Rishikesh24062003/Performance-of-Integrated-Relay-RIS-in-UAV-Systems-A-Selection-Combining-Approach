% clear all;
close all;
clc;

%% Parameters for Nakagami-m fading channel with N RIS elements
N = 16; % Number of reflecting elements (IRS)
M = 1e6; % Number of symbols (bits)
Average_SNR_dB = -50:5:10; % Average SNR in dB
Average_SNR = 10.^(Average_SNR_dB ./ 10); % Convert SNR to linear scale

m = 2; % Shape parameter for Nakagami-m distribution
omega = 1; % Scale parameter for Nakagami-m distribution

% Initialize arrays for SER simulation and analytical results
SER_sim = zeros(size(Average_SNR_dB));
P_e = zeros(size(Average_SNR_dB));

% Generate BPSK symbols
ip = rand(1, M) > 0.5; % Generate 0,1 with equal probability
s = 2 * ip - 1; % BPSK modulation: 0 -> -1, 1 -> 1 

% Simulation loop over different average SNR values
for jj = 1:length(Average_SNR_dB)
    % Calculate noise variance
    sigma = 1 / sqrt(2 * Average_SNR(jj));

    % Generate Nakagami-m fading coefficients for N elements
    h1 = sqrt(gamrnd(m, omega/m, N, M)); % |h|^2 is gamma distributed, |h| is sqrt(gamma)
    h2 = sqrt(gamrnd(m, omega/m, N, M)); % |g|^2 is gamma distributed, |g| is sqrt(gamma)
    
    % Calculate combined channel effect
    h = sum(h1 .* h2);

    % Add AWGN
    y = h .* s + sigma * (randn(1, M) + 1i * randn(1, M));

    % BPSK demodulation and Symbol Error Rate (SER) calculation
    ipHat = real(y) > 0;
    SER_sim(jj) = sum(ip ~= ipHat) / M;
end

% Analytical SER Calculation (BPSK) for IRS
lambda = sqrt((m * m) / (omega * omega));
E_vk1 = (gamma(m + 0.5) * gamma(m + 0.5)) * lambda^(-1) / (gamma(m) * gamma(m));
E_vk2 = (gamma(m + 1) * gamma(m + 1)) * lambda^(-2) / (gamma(m) * gamma(m));
ak = (N * (E_vk1)^2) / (E_vk2 - (E_vk1^2));
bk = (E_vk1) / (E_vk2 - (E_vk1^2));

p = 1;
q = 2; % for BPSK
z = 20; % the Gauss-Chebyshev integral parameter

% Analytical SER calculation loop
for j = 1:length(Average_SNR_dB)
    for i = 1:z
        x = cos((2 * i - 1) * pi / (2 * z));
        pe(i) = sqrt(1 - x^2) * pi / z * (1 / sqrt(2 * q * (log(2) - log(1 + x)))) * gammainc((bk / sqrt(Average_SNR(j)) * sqrt(2 * (log(2) - log(1 + x)) / q)), ak, 'lower');
    end
    P_e(j) = p * sqrt(q) / (2 * sqrt(2 * pi)) * sum(pe);
end

% Plotting SER results
figure;
semilogy(Average_SNR_dB, SER_sim, 'r-s', 'LineWidth', 1.5);
hold on;
semilogy(Average_SNR_dB, P_e, 'k--', 'LineWidth', 1.5);
xlabel('Average SNR (Es/N0) / dB');
ylabel('Symbol Error Rate (SER)');
title('SER vs Average SNR for Nakagami-m Fading Channel with IRS');
legend('Simulation', 'Analytical');
grid on;
axis([-50 10 1e-6 1])