clear;
close all;

%% Parameters
M = 1e6; % Number of symbols
Average_SNR_dB = -30:5:50; % Average SNR in dB
Average_SNR = 10.^(Average_SNR_dB ./ 10); % Convert SNR to linear scale

% Nakagami-m parameters
m = 2; % Shape parameter
Omega_SR= 1;
Omega_RD= 1; % Omega parameter for Nakagami-m distribution

% Modulation parameters (generalized)
modOrder = 2; % BPSK modulation order
s = 1; % General parameter for analytical expression
t = 2; % General parameter for analytical expression

% Free Space Path Loss parameters
d1 = 10; % Distance from source to relay in meters
d2 = 10; % Distance from relay to destination in meters
f = 1.6e9; % Frequency in Hz
c = 3e8; % Speed of light in m/s
lambda = c / f; % Wavelength

% Gains
G_s = 10^(10 / 10);
G_r = 10^(10 / 10);
G_d = 10^(10 / 10);

% Path Losses
PL1 = (G_s * G_r * lambda) / (4 * pi * d1^2);
PL2 = (G_d * G_r * lambda) / (4 * pi * d2^2);

% Generate symbols for any modulation order (assuming BPSK for now)
ip = randi([0, modOrder-1], 1, M); % Generate symbols
s = 2 * ip - 1; % BPSK modulation 0 -> -1, 1 -> 1 

% Initialize error count array
nErr = zeros(1, length(Average_SNR_dB));

for jj = 1:length(Average_SNR_dB)
    % Generate Nakagami-m fading coefficients
    h1 = sqrt(Omega_SR/2/m) * (randn(1, M) + 1j*randn(1, M)) .* sqrt(gamrnd(m, 1, [1, M]));
    h2 = sqrt(Omega_RD/2/m) * (randn(1, M) + 1j*randn(1, M)) .* sqrt(gamrnd(m, 1, [1, M]));

    % Apply path loss
    h1f = PL1 * h1;
    h2f = PL2 * h2;

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
Sim_Ser = nErr / M;

% SER Calculation (Analytical)

% Gauss-Laguerre quadrature parameters
p=1;
q=2;
% Compute zeros and weights for the Gauss-Laguerre quadrature
n=30;
alpha=-1/2;
[x,wl]=gauss_laguerre_weights(n, alpha);

for ii = 1:length(Average_SNR_dB) 

     summ=0;
   for t=1:n 
    th_sr = x(t)*Omega_SR;  % Normalized threshold for Nakagami-m CDF
    F_gamma_SR = 1-gammainc(m *(1/sqrt(PL1))* th_sr./Average_SNR(ii), m, 'lower');

    th_rd = x(t)*Omega_RD;  % Normalized threshold for Nakagami-m CDF
    F_gamma_RD = 1-gammainc(m *(1/sqrt(PL2))* th_rd./Average_SNR(ii), m, 'lower');

   summ=summ+((wl(t)*p.*sqrt(q))./(2*sqrt(2*pi))).*(1-(F_gamma_SR.*F_gamma_RD));
   end
   BER_any(ii)=summ;
end

%% Plotting SER results
figure;
semilogy(Average_SNR_dB, Sim_Ser, 'r-s', 'LineWidth', 1.5);
hold on;
semilogy(Average_SNR_dB, BER_any, 'k--', 'LineWidth', 1.5);
xlabel('Average SNR (dB)');
ylabel('Symbol Error Rate (SER)');
title(['SER vs Average SNR for Modulation Scheme in Nakagami-m Fading Channel']);
legend('Simulation','Analytical');
grid on;
% axis([-10 30 1e-6 1])
%% Function for Gauss-Laguerre quadrature
function [x, w] = gauss_laguerre_weights(n, alpha)
    i = 1:n;
    a = (2*i-1) + alpha;
    b = sqrt(i(1:n-1) .* ((1:n-1) + alpha));
    CM = diag(a) + diag(b,1) + diag(b,-1);
    [V, L] = eig(CM);
    [x, ind] = sort(diag(L));
    V = V(:, ind)';
    w = gamma(alpha+1) .* V(:,1).^2;
end 
