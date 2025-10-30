%% Parameters for Outage Probability Calculation
M = 100000;          % Number of symbols
Pt_dB = -20:5:50;      % Transmit power in dB
Pt = 10.^(Pt_dB/10); % Transmit power in linear scale
No = 1;              % Noise power
m = 2;               % Nakagami-m fading parameter
omega = 1;           % Omega parameter for Nakagami-m distribution

% Threshold for outage probability
SNR_th_dB = 10;      % Threshold SNR in dB
SNR_th = 10^(SNR_th_dB/10);

% Generate BPSK symbols
ip = rand(1, M) > 0.5; % Generate 0,1 with equal probability
s = 2 * ip - 1;        % BPSK modulation 0 -> -1, 1 -> 1

% Initialize outage probability array
OutageProb = zeros(size(Pt_dB));

%% Initialize 1D Network Parameters
NX = 200;               % Scenario X Length (1D space)
MAXS = 10;              % Maximum Speed
MAXT = 5;               % Maximum Moving Time
NoN = 2;               % Number of Nodes
simulationtime = 100;  % Number of Iterations

% Initialize matrices and arrays
graph = zeros(NX, 1);               % Initialize graph matrix for 1D
location = zeros(NoN, 3);           % Initialize location matrix
naibor_per_round = zeros(NoN, simulationtime); % Initialize neighbor per round matrix
naiborlocation = zeros(2, NoN);     % Initialize neighbor location matrix (1D)

% Source at the beginning
location(1,:) = [1, 0, 0];

% Destination at the end
location(NoN,:) = [NX, 0, 0];

% Randomly position UAV with relay and other nodes in between
for i = 2:NoN-1
    location(i,1) = randi(NX);
    if location(i,1) == 0
        location(i,1) = 1;
    end
    location(i,2) = 0; % Y-coordinate (not used in 1D)
    location(i,3) = 0; % Mobility parameter (not used in 1D)
    graph(location(i,1)) = graph(location(i,1)) + 1;
    naiborlocation(:, i) = [location(i,1); 0]; % Store neighbor location (1D)
end

%% 1D Random Waypoint (RWP) Mobility Model Simulation
s = 1; % Mobility model selector (1 for RWP)
for j = 1:simulationtime
    % Perform mobility update for each node
    for puse = 1:NoN
        if location(puse, 3) == 0
            % Node is active, add to work queue
            workque(puse) = puse; 
        else
            % Node is in pause, decrement pause time
            location(puse, 3) = location(puse, 3) - 1;  
        end
    end
    
    % Process nodes in the work queue
    for c = 1:NoN
        switch s
            case 1
                % RWP mobility model (simplified for 1D)
                if location(workque(c), 3) == 0
                    % Generate random movement within MAXS
                    moveDist = randi([-MAXS, MAXS]);
                    newLoc = location(workque(c), 1) + moveDist;
                    
                    % Ensure new location stays within bounds
                    if newLoc < 1
                        newLoc = 1;
                    elseif newLoc > NX
                        newLoc = NX;
                    end
                    
                    % Update graph and location
                    graph(location(workque(c), 1)) = graph(location(workque(c), 1)) - 1;
                    graph(newLoc) = graph(newLoc) + 1;
                    location(workque(c), 1) = newLoc;
                    
                    % Reset mobility parameter (not used in 1D)
                    location(workque(c), 3) = MAXT;
                end
        end
    end
    
    % Calculate outage probability for each Pt_dB
    for jj = 1:length(Pt_dB)
        % Generate Nakagami-m fading coefficients
        h1 = sqrt(omega/2/m) * (randn(1, M) + 1j*randn(1, M)) .* sqrt(gamrnd(m, 1, [1, M]));
        h2 = sqrt(omega/2/m) * (randn(1, M) + 1j*randn(1, M)) .* sqrt(gamrnd(m, 1, [1, M]));

        % Noise standard deviation
        sigma = 1 ./ (sqrt(2 * (Pt(jj) / No)));
        n1 = randn(1, M) + 1j * randn(1, M); % AWGN

        % Received signal at relay
        yR = h1 .* s + sigma .* n1;

        % Decode and forward
        sHat = real(yR ./ h1) > 0; % Equalize and make decision
        sHat = 2 * sHat - 1; % BPSK demodulation

        % Noise at destination
        n2 = randn(1, M) + 1j * randn(1, M);

        % Received signal at destination
        yD = h2 .* sHat + sigma .* n2;

        % Calculate instantaneous SNR
        SNR = Pt(jj) * min(abs(h1).^2, abs(h2).^2) / No;

        % Outage Probability
        OutageProb(jj) = OutageProb(jj) + mean(SNR < SNR_th);
    end
end

% Average outage probability over simulation time
OutageProb = OutageProb / simulationtime;

% Plotting Outage Probability results
figure;
semilogy(Pt_dB, OutageProb, 'r-s');
xlabel('Transmit Power (dB)');
ylabel('Outage Probability');
title('Outage Probability vs Transmit Power for Nakagami-m Fading Channel (1D RWP)');
grid on;