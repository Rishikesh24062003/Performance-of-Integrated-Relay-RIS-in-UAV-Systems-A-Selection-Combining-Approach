clear all;
clc;

% Parameters
m = 1;              % Nakagami-m shape factor (m >= 0.5)
omega = 1;          % Nakagami-m spread factor
num_samples = 1e6;  % Number of samples
snr_dB_range = 0:2:30; % SNR values from 0 to 30 dB in steps of 2 dB

% Free Space Path Loss parameters
d = 10; % Distance in meters
f = 1.6e9;  % Frequency in Hz (e.g., 2 GHz)
c = 3e8;  % Speed of light in m/s

% Calculate Free Space Path Loss (FSPL)
FSPL_dB = 20*log10(d) + 20*log10(f) + 20*log10(4*pi / c);
FSPL = 10^(FSPL_dB / 10); % Convert FSPL to linear scale

% Initialize array for SER
ser_values = zeros(size(snr_dB_range));

for i = 1:length(snr_dB_range)
    % Convert SNR to linear scale
    SNR_linear = 10^(snr_dB_range(i) / 10);

    % Generate Nakagami-m fading coefficients
    h = sqrt(gamrnd(m, omega/m, num_samples, 1));

    % Generate BPSK modulated transmitted signal
    tx_signal = randi([0 1], num_samples, 1) * 2 - 1; % BPSK: {0,1} -> {-1,1}

    % Apply Nakagami-m fading
    rx_signal = h .* tx_signal;

    % Generate complex Gaussian noise
    noise = (1/sqrt(2)) * (randn(num_samples, 1) + 1i*randn(num_samples, 1));

    % Add noise and FSPL to the received signal
    rx_signal_noisy = rx_signal + sqrt(FSPL / SNR_linear) * noise;

    % Demodulate received signal
    rx_demod = real(rx_signal_noisy) > 0;

    % Calculate SER
    ser_values(i) = mean(rx_demod ~= (tx_signal > 0));
end

% Plot SER vs SNR
figure;
semilogy(snr_dB_range, ser_values, 'r-o', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Symbol Error Rate (SER)');
title('Symbol Error Rate vs SNR for Nakagami-m Channel with FSPL');
legend(['Nakagami-m, m = ' num2str(m) ', \omega = ' num2str(omega)]);
