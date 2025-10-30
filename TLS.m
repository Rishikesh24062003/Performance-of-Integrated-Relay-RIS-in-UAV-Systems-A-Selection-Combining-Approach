clear all;
clc;

% Parameters
SNR_dB = 0:5:30;                  % SNR range in dB
N = 32;                          % Number of reflecting elements in LIS
M = 2;                            % M-QAM order
bpcu = log2(M);                   % bits per channel use (SE)
bits_total = bpcu * 8e5;          % Total number of bits (increase for higher SNR)
threshold_SNR_dB = 10;            % Outage threshold SNR in dB
threshold_SNR = 10^(threshold_SNR_dB / 10); % Outage threshold SNR in linear scale

% Preallocate arrays for results
outage_prob_intelligent = zeros(size(SNR_dB));
outage_prob_blind = zeros(size(SNR_dB));
SER_intelligent = zeros(size(SNR_dB));
SER_blind = zeros(size(SNR_dB));

for snr_idx = 1:length(SNR_dB)
    current_SNR_dB = SNR_dB(snr_idx);
    current_SNR = 10^(current_SNR_dB / 10);
    sigma = sqrt(1 / (2 * current_SNR)); % Noise variance
    bit_errors_intelligent = 0;
    bit_errors_blind = 0;
    outage_count_intelligent = 0;
    outage_count_blind = 0;

    % Generate random bit sequence
    bit_seq = randi([0 1], 1, bits_total);
    decoded_intelligent = zeros(1, bits_total);
    decoded_blind = zeros(1, bits_total);

    % Symbol generation
    ss = pskmod(0:M-1, M, 0, 'Gray'); % Gray-encoded M-PSK signal constellation
    ss = ss ./ sqrt(mean(abs(ss).^2)); % Normalization

    kk = 1;
    for iteration = 1:(bits_total / bpcu)
        bits = bit_seq(kk:kk + bpcu - 1);

        % Bin2Dec & Symbols
        dec = bi2de(bits, 'left-msb') + 1;
        x = ss(dec);

        % Channel generation
        h = (randn(N, 1) + 1i * randn(N, 1)) / sqrt(2);
        g = (randn(N, 1) + 1i * randn(N, 1)) / sqrt(2);

        % Intelligent Transmission
        Phi_intelligent = exp(-1i * angle(h .* g));
        B_intelligent = sum(h .* Phi_intelligent .* g);
        n = sigma * (randn + 1i * randn);
        r_intelligent = B_intelligent * x + n;
        SNR_intelligent = abs(B_intelligent)^2 / sigma^2;
        if SNR_intelligent < threshold_SNR
            outage_count_intelligent = outage_count_intelligent + 1;
        end
        [~, detected_intelligent] = min(arrayfun(@(k) norm(r_intelligent - B_intelligent * ss(k))^2, 1:M));
        decoded_intelligent(kk:kk + bpcu - 1) = de2bi(detected_intelligent - 1, bpcu, 'left-msb');

        % Blind Transmission
        %Phi_blind = ones(N, 1);
        % B_blind = sum(h .* Phi_blind .* g);
        % r_blind = B_blind * x + n;
        % SNR_blind = abs(B_blind)^2 / sigma^2;
        % if SNR_blind < threshold_SNR
        %     outage_count_blind = outage_count_blind + 1;
        % end
        % [~, detected_blind] = min(arrayfun(@(k) norm(r_blind - B_blind * ss(k))^2, 1:M));
        % decoded_blind(kk:kk + bpcu - 1) = de2bi(detected_blind - 1, bpcu, 'left-msb');
        % 
        % kk = kk + bpcu;
    end

    % Calculate bit errors
    bit_errors_intelligent = sum(bit_seq ~= decoded_intelligent);
    bit_errors_blind = sum(bit_seq ~= decoded_blind);

    % Calculate SER and outage probability
    SER_intelligent(snr_idx) = bit_errors_intelligent / bits_total;
    SER_blind(snr_idx) = bit_errors_blind / bits_total;
    outage_prob_intelligent(snr_idx) = outage_count_intelligent / (bits_total / bpcu);
    outage_prob_blind(snr_idx) = outage_count_blind / (bits_total / bpcu);
end

% Plot results
figure;
subplot(2, 1, 1);
semilogy(SNR_dB, outage_prob_intelligent, 'b-o', 'LineWidth', 1.5);
hold on;
semilogy(SNR_dB, outage_prob_blind, 'r-s', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Outage Probability');
legend('Intelligent Transmission', 'Blind Transmission');
title('Outage Probability vs SNR');

subplot(2, 1, 2);
semilogy(SNR_dB, SER_intelligent, 'b-o', 'LineWidth', 1.5);
hold on;
semilogy(SNR_dB, SER_blind, 'r-s', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Symbol Error Rate');
legend('Intelligent Transmission', 'Blind Transmission');
title('Symbol Error Rate vs SNR');
