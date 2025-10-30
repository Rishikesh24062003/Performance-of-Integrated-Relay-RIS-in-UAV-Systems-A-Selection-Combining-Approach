% Define Nakagami-m parameters
m = 1;        % Shape factor (m >= 0.5)
omega = 1;    % Spread factor

% Number of samples
num_samples = 1000000;

% Define SNR range in dB
snr_dB_range = 0:2:30; % SNR values from 0 to 30 dB in steps of 2 dB

% Initialize array for SER
ser_values = zeros(size(snr_dB_range));

% Loop over each SNR value to calculate SER
for i = 1:length(snr_dB_range)
    % Extract current SNR value in dB
    SNR_dB = snr_dB_range(i);
    
    % Convert SNR to linear scale
    SNR_linear = 10^(SNR_dB/10);
    
    % Generate Nakagami-m distributed channel coefficients
    h = sqrt(gamrnd(m, omega/m, num_samples, 1));
    
    % Generate BPSK modulated transmitted signal
    tx_signal = randi([0 1], num_samples, 1) * 2 - 1; % BPSK: {0,1} -> {-1,1}
    
    % Apply Nakagami-m fading
    rx_signal = h .* tx_signal;
    
    % Generate complex Gaussian noise
    noise = (1/sqrt(2)) * (randn(num_samples, 1) + 1i*randn(num_samples, 1));
    
    % Add noise to the received signal
    rx_signal_noisy = rx_signal + (1/sqrt(SNR_linear)) * noise;
    
    % Demodulate received signal
    rx_demod = real(rx_signal_noisy) > 0;
    
    % Calculate SER
    ser_values(i) = sum(rx_demod ~= (tx_signal > 0)) / num_samples;
end

% Plot SER vs SNR
figure;
semilogy(snr_dB_range, ser_values, 'r-o', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Symbol Error Rate (SER)');
title('Symbol Error Rate vs SNR for Nakagami-m Channel');
legend(['Nakagami-m, m = ' num2str(m) ', \omega = ' num2str(omega)]);
