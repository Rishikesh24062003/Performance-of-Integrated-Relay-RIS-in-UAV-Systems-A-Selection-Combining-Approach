% Clear workspace, close figures, and clear command window
clearvars;
close all;
clc;

% Parameters for simulation
M = 1e5;             % Number of symbols
Pt_dB = -20:2:30;    % Transmit power in dB
Pt = 10.^(Pt_dB/10); % Transmit power in linear scale
No = 1;              % Noise power
SNR_th_dB = 10;      % Threshold SNR in dB
SNR_th = 10^(SNR_th_dB/10); % Threshold SNR in linear scale

% Relay parameters
num_relays = 5;      % Number of relay nodes
relay_positions = rand(num_relays, 3) * 100;  % Random positions for relays

% Simulation: Relay Model
source_position = [0, 0, 0];        % Position of the source node
destination_position = [100, 100, 0]; % Position of the destination node

% Calculate distances
dist_source_relay = sqrt(sum((relay_positions - source_position).^2, 2));
dist_relay_destination = sqrt(sum((destination_position - relay_positions).^2, 2));

% Path loss function
PL_source_relay = path_loss(dist_source_relay);
PL_relay_destination = path_loss(dist_relay_destination);

% Initialize outage probability array
OutageProb = zeros(size(Pt_dB));

% Compute Outage Probability
for i = 1:length(Pt)
    % SNR at relays
    SNR_source_relay = Pt(i) ./ (No * PL_source_relay);
    SNR_relay_destination = Pt(i) ./ (No * PL_relay_destination);
    
    % Total SNR
    SNR_total = SNR_source_relay .* SNR_relay_destination ./ (SNR_source_relay + SNR_relay_destination + 1);
    
    % Outage Probability
    OutageProb(i) = mean(SNR_total < SNR_th);
end

% Plot: Outage Probability vs Transmit Power
figure;
semilogy(Pt_dB, OutageProb, 'r-s', 'LineWidth', 1.5);
xlabel('Transmit Power (dB)');
ylabel('Outage Probability');
title('Outage Probability vs Transmit Power');
grid on;

% Path loss function
function PL = path_loss(distances)
    % Path loss parameters
    n = 3;
    B = [735, -1190, 455] / 72;
    beta = [2, 4, 6];
    D = max(distances);

    % Calculate path loss
    PL = zeros(size(distances));
    for j = 1:n
        PL = PL + (B(j) / D^(beta(j) + 1)) .* distances.^beta(j);
    end
end
