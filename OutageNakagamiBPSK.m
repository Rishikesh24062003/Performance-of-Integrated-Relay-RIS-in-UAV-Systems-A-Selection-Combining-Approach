% Define Nakagami-m parameters
m = 1; % shape factor (m >= 0.5)
omega = 1; % spread factor

% Number of samples
num_samples = 1000000;

% Define SNR range in dB
snr_dB_range = 0:5:30; % SNR values from 0 to 30 dB in steps of 2 dB
snr_dB_range_linear = 10.^(snr_dB_range./10);

% Define outage SNR threshold in dB
snr_threshold_dB = 5; % Outage threshold in dB
snr_threshold_linear = 10.^(snr_threshold_dB./10);


% Initialize array for outage probability
%outage_probabilities = zeros(size(snr_dB_range));

% Loop over each SNR value to calculate outage probability
for i = 1:length(snr_dB_range)

    SNR_dB = snr_dB_range(i);
    
    % Generate Nakagami-m distributed channel coefficients
    h = sqrt(gamrnd(m, omega/m, 1, num_samples));
    
    % Calculate received SNR
    instt_snr = (abs(h).^2) * snr_dB_range_linear(i);
    count = 0;
    for j =1:num_samples
        if instt_snr(j) < snr_threshold_linear
            count = count +1;
        end
            
    end
    % Calculate outage probability
    nErr = count./num_samples
    outage_probabilities(i) = nErr;
end

% Plot Outage Probability vs SNR
figure;
semilogy(snr_dB_range, outage_probabilities, 'b-o', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Outage Probabilities');
title('Outage Probability vs SNR for Nakagami-m Channel');
legend(['Nakagami-m, m = ' num2str(m) ', \omega = ' num2str(omega)]);
