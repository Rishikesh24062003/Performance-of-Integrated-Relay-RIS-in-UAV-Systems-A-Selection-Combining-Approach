clear all;
clc;

% Parameters
m = 1; % Shape factor (m >= 0.5)
omega = 1; % Spread factor
num_samples = 1e6; % Number of samples
snr_dB_range = -10:5:30; % SNR values from 0 to 30 dB in steps of 5 dB
snr_threshold_dB = -50; % Outage threshold in dB

% Convert SNR values to linear scale
snr_linear_range = 10.^(snr_dB_range ./ 10);
snr_threshold_linear = 10^(snr_threshold_dB / 10);

% Free Space Path Loss parameters
d = 10; % Distance in meters
f = 1.6e9; % Frequency in Hz (e.g., 2 GHz)
c = 3e8; % Speed of light in m/s

% Calculate Free Space Path Loss (FSPL)
FSPL_dB = 20*log10(d) + 20*log10(f) + 20*log10(4 * pi / c);
FSPL = 10^(FSPL_dB / 10); % Convert FSPL to linear scale

% Initialize array for outage probabilities
outage_probabilities = zeros(size(snr_dB_range));

% Loop over each SNR value to calculate outage probability
for i = 1:length(snr_dB_range)
    % Generate Nakagami-m distributed channel coefficients
    h = sqrt(gamrnd(m, omega/m, 1, num_samples));
    
    % Calculate received SNR with FSPL
    instt_snr = (abs(h).^2) * snr_linear_range(i) / FSPL;
    
    % Calculate outage probability
    outage_probabilities(i) = mean(instt_snr < snr_threshold_linear);
end

% Plot Outage Probability vs SNR
figure;
semilogy(snr_dB_range, outage_probabilities, 'b-o', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Outage Probability');
title('Outage Probability vs SNR for Nakagami-m Channel with FSPL');
legend(['Nakagami-m, m = ' num2str(m) ', \omega = ' num2str(omega)]);
