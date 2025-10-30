% Clear workspace, close figures, and clear command window
clearvars;
close all;
clc;

% Parameters for simulation
M = 1e5;             % Number of symbols
Pt_dB = -20:2:30;    % Transmit power in dB
Pt = 10.^(Pt_dB/10); % Transmit power in linear scale
No = 1;              % Noise power
m = 2;               % Nakagami-m fading parameter
omega = 1;           % Omega parameter for Nakagami-m distribution

% Relay parameters
num_relays = 5;      % Number of relay nodes
relay_positions = rand(num_relays, 3) * 100;  % Random positions for relays

% Simulation: Relay Model
source_position = [0, 0, 0];       % Position of the source node
destination_position = [100, 100, 0];  % Position of the destination node

% Calculate distances
dist_source_relay = sqrt(sum((relay_positions - source_position).^2, 2));
dist_relay_destination = sqrt(sum((destination_position - relay_positions).^2, 2));

% Calculate path losses
PL_source_relay = path_loss(dist_source_relay);
PL_relay_destination = path_loss(dist_relay_destination);

% Calculate SNRs at relays
SNR_source_relay = repmat(Pt', 1, num_relays) ./ (No * PL_source_relay');
SNR_relay_destination = repmat(Pt', 1, num_relays) ./ (No * PL_relay_destination');

% Calculate total SNR through relays using harmonic mean
SNR_total = 2 ./ (1 ./ SNR_source_relay + 1 ./ SNR_relay_destination);

% Compute Symbol Error Rate (SER) using Nakagami-m fading model
SER = zeros(size(Pt));
for i = 1:length(Pt)
    gamma_m = SNR_total(i, :);
    SER(i) = mean(1 - gammainc(m, gamma_m / omega));
end

% Plot: Symbol Error Rate (SER) vs Transmit Power
figure;
semilogy(Pt_dB, SER, 'b-o', 'LineWidth', 1.5);
xlabel('Transmit Power (dB)');
ylabel('Symbol Error Rate (SER)');
title('Symbol Error Rate (SER) vs Transmit Power');
grid on;

% Function for path loss calculation
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
