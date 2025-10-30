clear all;
close all;
clc;

%% Parameters
N = 16; % Number of reflecting elements (IRS)
M = 10^6; % Number of symbols (bits)
Average_SNR_dB = -50:5:10; % Average SNR in dB
Average_SNR = 10.^(Average_SNR_dB ./ 10); % Convert SNR to linear scale

m = 2; % Shape parameter for Nakagami-m distribution
omega = 1; % Scale parameter for Nakagami-m distribution
Omega_SR = 2;
Omega_RD = 2;

d1 = 5; % Distance from source to relay in meters
d2 = 5; % Distance from relay to destination in meters
f = 1.6e9; % Frequency in Hz
c = 3e8; % Speed of light in m/s
lambda = c / f; % Wavelength

G_s = 10^(20 / 10);
G_r = 10^(20 / 10);
G_d = 10^(20 / 10);

PL1 = (G_s * G_r * lambda) / (4 * pi * d1^2);
PL2 = (G_d * G_r * lambda) / (4 * pi * d2^2);

%% Simulation: IRS
% Initialize arrays for SER simulation and analytical results
SER_sim_IRS = zeros(size(Average_SNR_dB));
P_e_IRS = zeros(size(Average_SNR_dB));

% Generate BPSK symbols
ip = rand(1, M) > 0.5; % Generate 0,1 with equal probability
s = 2 * ip - 1; % BPSK modulation: 0 -> -1, 1 -> 1 

% Simulation loop over different average SNR values
for jj = 1:length(Average_SNR_dB)
    % Calculate noise variance
    sigma = 1 / sqrt(2 * Average_SNR(jj));

    % Simulation of Nakagami-m fading channel coefficients
    hi = sqrt(gamrnd(m, omega/m, N, M)); % Nakagami-m fading channel coefficients (hi)
    gi = sqrt(gamrnd(m, omega/m, N, M)); % Nakagami-m fading channel coefficients (gi)
    h1 =  hi .* gi; % Nakagami-m fading channel (h1)

    h = sum(h1); % Sum of channel effects over N elements

    % Add AWGN
    y = h .* s + sigma * (randn(1, M) + 1i * randn(1, M));

    % BPSK demodulation and Symbol Error Rate (SER) calculation
    ipHat = real(y) > 0;
    SER_sim_IRS(jj) = sum(ip ~= ipHat) / M;
end

%% Analytical SER: IRS
lambda = sqrt((m * m) / (omega * omega));
E_vk1 = (gamma(m + 0.5) * gamma(m + 0.5)) * lambda^(-1) / (gamma(m) * gamma(m));
E_vk2 = (gamma(m + 1) * gamma(m + 1)) * lambda^(-2) / (gamma(m) * gamma(m));
ak = (N * (E_vk1)^2) / (E_vk2 - (E_vk1^2));
bk = (E_vk1) / (E_vk2 - (E_vk1^2));

p = 1;
q = 2; % for BPSK
z = 20; % the Gauss-Chebyshev integral parameter

for j = 1:length(Average_SNR_dB)
    pe = zeros(1, z);
    for i = 1:z
        x = cos((2 * i - 1) * pi / (2 * z));
        pe(i) = sqrt(1 - x^2) * pi / z * (1 / sqrt(2 * q * (log(2) - log(1 + x)))) * ...
                gammainc((bk / sqrt(Average_SNR(j)) * sqrt(2 * (log(2) - log(1 + x)) / q)), ak, 'lower');
    end
    P_e_IRS(j) = p * sqrt(q) / (2 * sqrt(2 * pi)) * sum(pe);
end

%% Simulation: Relay
% Initialize error count array
nErr = zeros(1, length(Average_SNR_dB));

for jj = 1:length(Average_SNR_dB)
    % Generate Nakagami-m fading coefficients
    h1 = sqrt(Omega_SR/2/m) * (randn(1, M) + 1j*randn(1, M)) .* sqrt(gamrnd(m, 1, [1, M]));
    h2 = sqrt(Omega_RD/2/m) * (randn(1, M) + 1j*randn(1, M)) .* sqrt(gamrnd(m, 1, [1, M]));

    % Apply path loss
    h1f =  PL1 * h1;
    h2f =  PL2 * h2;

    % Noise standard deviation
    sigma = sqrt(1/ (2 * Average_SNR(jj)));

    % Generate noise
    n1 = sigma * (randn(1, M) + 1j * randn(1, M)); % AWGN

    % Received signal at relay and decode
    yR = h1f .* s + n1;
    sHat = real(yR ./ h1f) > 0; % Equalize and make decision
    sHat = 2 * sHat - 1; % BPSK demodulation

    % Generate noise at destination
    n2 = sigma * (randn(1, M) + 1j * randn(1, M));
    
    % Received signal at destination
    yD = h2f .* sHat + n2;
    ipHat = real(yD ./ h2f) > 0; % Equalization and decision

    % Count errors
    nErr(jj) = sum(ip ~= ipHat);
end

% Calculate SER
Sim_Ser_Relay = nErr / M;

%% Analytical SER: Relay
n = 30;
alpha = -0.5;
[x, wl] = gauss_laguerre_weights(n, alpha);

BER_any_Relay = zeros(size(Average_SNR_dB));
for ii = 1:length(Average_SNR_dB)
    summ = 0;
    for t = 1:n
        th_sr = x(t) * Omega_SR;
        F_gamma_SR = 1 - gammainc(m * (1 / sqrt(PL1)) * th_sr / Average_SNR(ii), m, 'lower');

        th_rd = x(t) * Omega_RD;
        F_gamma_RD = 1 - gammainc(m * (1 / sqrt(PL2)) * th_rd / Average_SNR(ii), m, 'lower');

        summ = summ + ((wl(t) * p .* sqrt(q)) / (2 * sqrt(2 * pi))) * (1 - (F_gamma_SR .* F_gamma_RD));
    end
    BER_any_Relay(ii) = summ;
end

%% Selection Combining
SER_SC_sim = SER_sim_IRS .* Sim_Ser_Relay;
SER_SC_analytical = P_e_IRS .* BER_any_Relay;

%% Plotting the results for Selection Combining
figure;
semilogy(Average_SNR_dB, SER_SC_sim, 'bo-', 'LineWidth', 2); % Blue circles for simulation
hold on;
semilogy(Average_SNR_dB, SER_SC_analytical, 'r*-', 'LineWidth', 2); % Red stars for analytical
grid on;
legend('Selection Combining Simulation', 'Selection Combining Analytical');
xlabel('Average SNR (dB)');
ylabel('Symbol Error Rate (SER)');
title('SER vs. Average SNR for Selection Combining');

%% Function for Gauss-Laguerre quadrature
function [x, w] = gauss_laguerre_weights(n, alpha)
    i = 1:n;
    a = (2*i-1) + alpha;
    b = sqrt(i(1:n-1) .* ((1:n-1) + alpha));
    CM = diag(a) + diag(b, 1) + diag(b, -1);
    [V, L] = eig(CM);
    [x, ind] = sort(diag(L));
    V = V(:, ind)';
    w = gamma(alpha + 1) .* V(:, 1).^2;
end
