% Parameters
M = 10000000; % Number of symbols
Pt_dB = -20:5:50; % Transmit power in dB
Pt = 10.^(Pt_dB/10); % Transmit power in linear scale
No = 1; % Noise power
m = 2; % Nakagami-m fading parameter
omega = 1; % Omega parameter for Nakagami-m distribution
SNR_th_dB = 10; % Threshold SNR in dB
SNR_th = 10^(SNR_th_dB/10); % Threshold SNR in linear scale
p_s = 0.5; % Probability that a node is static
p_p = 0.3; % Example value for p_p
f_X = @(x) p_s * exp(-x) + (1 - p_s) * (p_p + (1 - p_p)); % Node distribution function

% Generate BPSK symbols
s = 2 * (rand(1, M) > 0.5) - 1;

% Generate Nakagami-m fading coefficients
h1i = sqrt(omega/2/m) * (randn(1, M) + 1j*randn(1, M)) .* sqrt(gamrnd(m, 1, [1, M]));
h2i = sqrt(omega/2/m) * (randn(1, M) + 1j*randn(1, M)) .* sqrt(gamrnd(m, 1, [1, M]));

% Generate random CDF values
a = 4;
n = M; % Set n to M to match the size
q = 1/8 - rand(1, n) / 4;
d = sqrt(1/64 - q.^2);
z = (q + 1j * d).^(1/3);
xrand = a * (1/2 - real(z) + imag(z) * sqrt(3));

% Initialize outage probability array
OutageProb = zeros(1, length(Pt_dB));

for jj = 1:length(Pt)
    h1 = h1i .* xrand;
    h2 = h2i .* xrand;

    % Noise
    sigma = 1 ./ sqrt(2 * Pt(jj) / No);
    n1 = randn(1, M) + 1j * randn(1, M);
    n2 = randn(1, M) + 1j * randn(1, M);

    % Positions according to node distribution
    x_sr = f_X(rand(1, M));
    x_rd = f_X(rand(1, M));

    % Received signal at relay and destination
    yR = h1 .* s + sigma .* n1;
    sHat = real(yR ./ h1) > 0;
    sHat = 2 * sHat - 1;
    yD = h2 .* sHat + sigma .* n2;

    % Instantaneous SNR for both links
    SNR1 = Pt(jj) * abs(h1).^2 ./ (No * x_sr.^2);
    SNR2 = Pt(jj) * abs(h2).^2 ./ (No * x_rd.^2);
    SNR_min = min(SNR1, SNR2);

    % Outage Probability
    OutageProb(jj) = mean(SNR_min < SNR_th);
end

% Plot results
figure;
semilogy(Pt_dB, OutageProb, 'r-s');
xlabel('Transmit Power (dB)');
ylabel('Outage Probability');
title('Outage Probability vs Transmit Power for Nakagami-m Fading Channel with Node Distribution');
grid on;
