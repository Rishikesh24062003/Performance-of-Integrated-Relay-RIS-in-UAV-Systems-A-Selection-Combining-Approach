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
SNR_th_dB = 10;      % Threshold SNR in dB
SNR_th = 10^(SNR_th_dB/10); % Threshold SNR in linear scale

% Random Waypoint Model parameters (example values)
area_size = 100;     % Size of the cubic area (in meters)
num_nodes = 10;      % Number of nodes
max_speed = 2;       % Maximum speed of the nodes (m/s)
pause_time = 1;      % Pause time at each waypoint (s)
sim_time = 100;      % Total simulation time (s)
delta_t = 0.1;       % Simulation time step (seconds)

% Initialize array to store Symbol Error Rate (SER) values
SymbolErrorRate = zeros(size(Pt));

% Simulation loop for each transmit power level
for k = 1:length(Pt)
    % Reset positions and velocities for each iteration
    positions = rand(num_nodes, 3) * area_size;
    destinations = rand(num_nodes, 3) * area_size;
    velocities = zeros(num_nodes, 3);
    pause_timer = zeros(num_nodes, 1);
    
    % Initialize array to store positions over time
    num_steps = sim_time / delta_t; % Number of simulation steps
    node_positions = zeros(num_nodes, 3, num_steps); % Preallocate position matrix
    
    % Simulation loop
    for t = 1:num_steps
        % Update positions and velocities for each node
        for i = 1:num_nodes
            if pause_timer(i) > 0
                % Node is currently pausing
                pause_timer(i) = pause_timer(i) - delta_t; % Decrement pause timer
            else
                % Update position based on velocity
                positions(i, :) = positions(i, :) + velocities(i, :) * delta_t;
                
                % Check if node has reached its destination
                if norm(positions(i, :) - destinations(i, :)) < 1
                    % Node reaches destination, choose new destination and speed
                    destinations(i, :) = rand(1, 3) * area_size; % New random destination
                    speed = rand() * max_speed; % New random speed
                    direction = (destinations(i, :) - positions(i, :)) / norm(destinations(i, :) - positions(i, :));
                    velocities(i, :) = direction * speed; % Update velocity
                    pause_timer(i) = pause_time; % Set pause time
                end
            end
        end
        
        % Store positions at each time step
        node_positions(:, :, t) = positions;
    end
    
    % Calculate distances between nodes and an arbitrary reference point (e.g., center of the area)
    reference_point = area_size / 2 * ones(1, 3);
    distances = sqrt(sum((reshape(node_positions, [], 3) - reference_point).^2, 2));
    
    % Function for path loss calculation (f_R(r))
    n = 3;
    B = [735, -1190, 455] / 72;
    beta = [2, 4, 6];
    D = max(distances);
    
    % Calculate path loss (f_R(r)) for each distance
    PL_i = f_R(distances, n, B, beta, D);
    
    % Compute instantaneous SNR for each Pt
    inst_SNR = Pt(k) * omega ./ (No * PL_i); % Adjusted for path loss PL_i
    
    % Ensure inst_SNR is real and non-negative
    inst_SNR = max(inst_SNR, eps); % Avoid inst_SNR being zero
    
    % Compute Symbol Error Rate (SER) using Q-function approximation for Nakagami-m fading
    SymbolErrorRate(k) = mean(1 - gammainc(m, inst_SNR / m));
end

% Plot: Symbol Error Rate (SER) vs Transmit Power (Pt_dB)
figure;
semilogy(Pt_dB, SymbolErrorRate, 'b-o', 'LineWidth', 1.5);
xlabel('Transmit Power (dB)');
ylabel('Symbol Error Rate (SER)');
title('Symbol Error Rate (SER) vs Transmit Power');
grid on;

% Function for path loss calculation (f_R(r))
function PL_i = f_R(distances, n, B, beta, D)
    if size(distances, 1) ~= 1
        distances = distances'; % Ensure distances is a row vector
    end
    PL_i = zeros(size(distances));
    for i = 1:n
        PL_i = PL_i + (B(i) / D^(beta(i) + 1)) .* distances.^beta(i);
    end
end
