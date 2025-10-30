%% Clear and close all
clear;
close all;
clc;

%% Parameters
M = 10^6; % Number of symbols
Average_SNR_dB = -30:5:25; % Average SNR in dB
Average_SNR = 10.^(Average_SNR_dB ./ 10); % Convert SNR to linear scale

% Nakagami-m parameters
m = 2; % Shape parameter
Omega = 1; % Omega parameter for Nakagami-m distribution

% Modulation parameters for BPSK
s = 1;
t = 2;

% Gauss-Laguerre quadrature parameters
n = 20; % Number of quadrature points
alpha = -0.5;

% Compute zeros and weights for the Gauss-Laguerre quadrature
[phi, W] = gauss_laguerre_weights(n, alpha);

%% ASER Calculation (Analytical)
Pser_analytical = zeros(size(Average_SNR));

for k = 1:length(Average_SNR)
    Omega = Average_SNR(k);
    % Sum contributions from each quadrature point
    sum_terms = 0;
    for l = 1:n
        F_gamma = 1 - exp(-m * phi(l) / Omega) * sum(arrayfun(@(i) (m * phi(l) / Omega)^i / factorial(i), 0:m-1));
        sum_terms = sum_terms + W(l) * F_gamma;
    end
    Pser_analytical(k) = (s * sqrt(t)) / (2 * sqrt(2 * pi)) * sum_terms;
end

%% Simulation Parameters
Pt_dB = Average_SNR_dB; % Use the same SNR range for simulation
Pt = Average_SNR; % Transmit power in linear scale
No = 1; % Noise power

% Free Space Path Loss parameters
d1 = 10; % Distance from source to relay in meters
d2 = 10; % Distance from relay to destination in meters
f = 1.6e9; % Frequency in Hz
c = 3e8; % Speed of light in m/s
lambda = c / f; % Wavelength

% Gains
G_s = 1; G_r = 1; G_d = 1;

% Path Losses
PL1 = (G_s * G_r * lambda) / (4 * pi * d1^2);
PL2 = (G_d * G_r * lambda) / (4 * pi * d2^2);

% Generate BPSK symbols
ip = rand(1, M) > 0.5; % Generate 0,1 with equal probability
s = 2 * ip - 1; % BPSK modulation 0 -> -1, 1 -> 1 

% Initialize error count array
nErr = zeros(1, length(Pt_dB));

for jj = 1:length(Pt)
    % Generate Nakagami-m fading coefficients
    h1 = sqrt(Omega/2/m) * (randn(1, M) + 1j*randn(1, M)) .* sqrt(gamrnd(m, 1, [1, M]));
    h2 = sqrt(Omega/2/m) * (randn(1, M) + 1j*randn(1, M)) .* sqrt(gamrnd(m, 1, [1, M]));

    % Apply path loss
    h1f = PL1 * h1;
    h2f = PL2 * h2;

    % Noise standard deviation
    sigma = sqrt(No / (2 * Pt(jj)));

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
simSer = nErr / M;

%% Plotting ASER results
figure;
semilogy(Average_SNR_dB, Pser_analytical, 'b-o', 'LineWidth', 1.5);
hold on;
% semilogy(Pt_dB, simSer, 'r-s', 'LineWidth', 1.5);
xlabel('Average SNR (dB)');
ylabel('Symbol Error Rate (SER)');
title(['SER vs Average SNR for BPSK in Nakagami-m Fading Channel']);
legend( 'Simulation');
grid on;

%% Functions for Gauss-Laguerre quadrature
function [x, w] = gauss_laguerre_weights(n, alpha)
    % This function determines the abscisas (x) and weights (w) for the
    % Gauss-Laguerre quadrature of order n>1, on the interval [0, +infinity].
    % Unlike the function 'GaussLaguerre', this function is valid for
    % n>=34. This is due to the fact that the companion matrix (of the n'th
    % degree Laguerre polynomial) is now constructed as a symmetrical
    % matrix, guaranteeing that all the eigenvalues (roots) will be real.

    % © Geert Van Damme
    % geert@vandamme-iliano.be
    % February 21, 2010

    % Building the companion matrix CM
    % CM is such that det(xI-CM)=L_n(x), with L_n the Laguerree polynomial
    % under consideration. Moreover, CM will be constructed in such a way
    % that it is symmetrical.
    i = 1:n;
    a = (2*i-1) + alpha;
    b = sqrt(i(1:n-1) .* ((1:n-1) + alpha));
    CM = diag(a) + diag(b,1) + diag(b,-1);

    % Determining the abscissas (x) and weights (w)
    % - since det(xI-CM)=L_n(x), the abscissas are the roots of the
    %   characteristic polynomial, i.d. the eigenvalues of CM;
    % - the weights can be derived from the corresponding eigenvectors.
    [V, L] = eig(CM);
    [x, ind] = sort(diag(L));
    V = V(:, ind)';
    w = gamma(alpha+1) .* V(:,1).^2;
end
