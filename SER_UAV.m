clear all;
clc;

% Parameters
M = 10^6; % Number of symbols
Pt_dB = 0:2:30; % Transmit power in dB
Pt = 10.^(Pt_dB/10); % Transmit power in linear scale
No = 1; % Noise power
m = 2; % Nakagami-m fading parameter
omega = 1; % Omega parameter for Nakagami-m distribution
N = 4; % Number of reflecting elements in IRS

% Generate BPSK symbols
ip = rand(1, M) > 0.5; % Generate 0,1 with equal probability
s = 2 * ip - 1; % BPSK modulation 0 -> -1, 1 -> 1 

% Initialize the SER array
SER = zeros(1, length(Pt_dB));

for jj = 1:length(Pt)
    % Generate Nakagami-m fading coefficients for relay
    % h_relay = sqrt(omega/2/m) * (randn(1, M) + 1j*randn(1, M)) .* sqrt(gamrnd(m, 1, [1, M]));
    h1 = sqrt(omega/2/m) * (randn(1, M) + 1j*randn(1, M)) .* sqrt(gamrnd(m, 1, [1, M]));
    h2 = sqrt(omega/2/m) * (randn(1, M) + 1j*randn(1, M)) .* sqrt(gamrnd(m, 1, [1, M]));
    h_relay = min(h1,h2);
    
    % Generate Nakagami-m fading coefficients for IRS with N elements
    hi = sqrt(gamrnd(m, omega/m, [N, M]));  % |h|^2 is gamma distributed, |h| is sqrt(gamma)
    gi = sqrt(gamrnd(m, omega/m, [N, M]));  % |g|^2 is gamma distributed, |g| is sqrt(gamma)
    h_IRS = sum(hi .* gi); % Combined IRS channel

    % Select the better channel for communication
    h_final = max(abs(h_relay), abs(h_IRS));

    % Noise standard deviation
    sigma = sqrt(No / (2 * Pt(jj)));

    % Generate noise
    n = sigma * (randn(1, M) + 1j * randn(1, M)); % AWGN

    % Received signal
    r = h_final .* s + n;

    % BPSK demodulation
    sHat = real(r ./ h_final) > 0; % Equalize and make decision
    sHat = 2 * sHat - 1; % BPSK demodulation

    % Calculate SER
    SER(jj) = mean(sHat ~= s);
end

% Plotting SER results
figure;
semilogy(Pt_dB, SER, 'b-o', 'LineWidth', 1.5);
xlabel('Transmit Power (dB)');
ylabel('Symbol Error Rate');
title('Symbol Error Rate vs Transmit Power for Nakagami-m Fading Channel with UAV (Relay and IRS)');
grid on;
