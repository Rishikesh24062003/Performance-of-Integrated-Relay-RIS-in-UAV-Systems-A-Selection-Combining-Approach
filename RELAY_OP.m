% Parameters (set as per your requirements)
M = 1e7;         % Number of symbols
Pt_dB = -10:5:40;  % Transmit power in dB
Pt = 10.^(Pt_dB/10); % Transmit power in linear scale
No = 1;          % Noise power
m = 2;           % Nakagami-m fading parameter
omega = 1;       % Omega parameter for Nakagami-m distribution

% Free Space Path Loss parameters
d1 = 10;         % Distance from source to relay in meters
d2 = 10;         % Distance from relay to destination in meters
f = 1.6e9;       % Frequency in Hz (1.6 GHz)
c = 3e8;         % Speed of light in m/s
lambda = c / f;  % Wavelength

% Gains (assuming unity gains for simplicity)
G_s = 1; G_r = 1; G_d = 1;

% Path Losses (Free Space Path Loss)
PL1 = (G_s * G_r * lambda) / (4 * pi * d1^2);
PL2 = (G_d * G_r * lambda) / (4 * pi * d2^2);

% Threshold for outage probability
SNR_th_dB = -10;  % Threshold SNR in dB
SNR_th = 10^(SNR_th_dB/10);

% Generate BPSK symbols
ip = rand(1, M) > 0.5; % Generate 0,1 with equal probability
s = 2 * ip - 1; % BPSK modulation 0 -> -1, 1 -> 1 

% Initialize the outage probability array
OutageProb = zeros(1, length(Pt_dB));

% Simulation loop over different transmit powers
for jj = 1:length(Pt)
    % Generate Nakagami-m fading coefficients
    h1 = sqrt(gamrnd(m, omega/m, 1, M)); % Nakagami-m fading coefficients for h1
    h2 = sqrt(gamrnd(m, omega/m, 1, M)); % Nakagami-m fading coefficients for h2

    % Apply path loss
    h1f = sqrt(PL1) * h1;
    h2f = sqrt(PL2) * h2;

    % Noise standard deviation adjusted for FSPL
    sigma = sqrt(No / (2 * Pt(jj)));

    % Generate noise
    n1 = sigma * (randn(1, M) + 1j * randn(1, M)); % AWGN

    % Received signal at relay
    yR = h1f .* s + n1;

    % Decode and forward
    sHat = real(yR ./ h1f) > 0; % Equalize and make decision
    sHat = 2 * sHat - 1; % BPSK demodulation

    % Noise at destination
    n2 = sigma * (randn(1, M) + 1j * randn(1, M));

    % Received signal at destination
    yD = h2f .* sHat + n2;

    % Calculate instantaneous SNR
    SNR1 = Pt(jj) * abs(h1f).^2 / No;
    SNR2 = Pt(jj) * abs(h2f).^2 / No;
    SNR = min(SNR1, SNR2); % Effective SNR is the minimum of the two links

    % Calculate outage probability
    OutageProb(jj) = mean(SNR < SNR_th);
end

% Plotting Outage Probability results
figure;
semilogy(Pt_dB, OutageProb, 'r-s', 'LineWidth', 1.5);
xlabel('Transmit Power (dB)');
ylabel('Outage Probability');
title('Outage Probability vs Transmit Power for Nakagami-m Fading Channel with FSPL');
grid on;
