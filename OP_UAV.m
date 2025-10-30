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
N = 32;              % Number of reflecting elements in IRS
SNR_th_dB = 10;      % Threshold SNR in dB
SNR_th = 10^(SNR_th_dB/10); % Threshold SNR in linear scale

% Random Waypoint Model parameters (example values)
area_size = 100;     % Size of the cubic area (in meters)
num_nodes = 10;      % Number of nodes
max_speed = 2;       % Maximum speed of the nodes (m/s)
pause_time = 1;      % Pause time at each waypoint (s)
sim_time = 100;      % Total simulation time (s)

% Simulation: Random Waypoint Model
positions = rand(num_nodes, 3) * area_size; % Random initial positions in the area
destinations = rand(num_nodes, 3) * area_size; % Random destinations for each node
velocities = zeros(num_nodes, 3); % Initial velocities (zero for now)
pause_timer = zeros(num_nodes, 1); % Timer for pause time
delta_t = 1; % Simulation time step (seconds)

% Initialize array to store positions over time
num_steps = sim_time / delta_t; % Number of simulation steps
node_positions = zeros(num_nodes * num_steps, 3); % Preallocate position matrix

% Simulation loop
current_time = 0; % Initialize current simulation time
step = 1; % Step counter
while current_time < sim_time
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
    node_positions((step - 1) * num_nodes + (1:num_nodes), :) = positions;
    
    % Update simulation time and step counter
    current_time = current_time + delta_t;
    step = step + 1;
end

% Calculate distances between nodes and an arbitrary reference point (e.g., center of the area)
reference_point = area_size / 2 * ones(1, 3);
distances = sqrt(sum((node_positions - reference_point).^2, 2));

% Function for path loss calculation (f_R(r))
n = 3;
B = [735, -1190, 455] / 72;
beta = [2, 4, 6];
D = max(distances);

% Calculate path loss (f_R(r)) for each distance
PL_i = f_R(distances, n, B, beta, D);

% Example plot: Outage Probability vs Transmit Power
OutageProb = rand(size(Pt_dB)); % Example data

% Plotting Outage Probability results
figure;
semilogy(Pt_dB, OutageProb, 'r-s', 'LineWidth', 1.5);
xlabel('Transmit Power (dB)');
ylabel('Outage Probability');
title('Outage Probability vs Transmit Power');
grid on;

% Display the graph immediately
drawnow;

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
