function [SOC_Estimated, Vt_Estimated, Vt_Error] = i_write4_aukf_2(Current, Vt_Actual, Temperature, RecordingTime)
% ASRUKF_SOC_Estimation - Improved Adaptive Square Root Unscented Kalman Filter for SOC Estimation
% Fixed issues with SOC tracking and initial state estimation

%% Load Parameters
load 'BatteryModel_AUKF.mat';    % loads `param` table among others
load 'SOC_OCV_AUKF.mat';         % loads `SOC_OCV_AUKF` table

%% Battery Constants
QN = 100 * 3600;    % amp-seconds (100 Ah battery)
eta = 1.0;          % Coulombic efficiency

% Calculate sampling time
%if numel(RecordingTime) > 1
   % deltaT = median(diff(RecordingTime));
%else
    deltaT = 30;
%end
fprintf('Calculated sampling time: %.4f seconds\n', deltaT);

%% Convert SOC_R strings (e.g. '100-90') into numeric midpoints [0,1]
numSOC = zeros(height(param), 1);
for i = 1:height(param)
    parts = regexp(param.SOC_R{i}, '-', 'split');
    upper = str2double(parts{1});
    lower = str2double(parts{2});
    numSOC(i) = (upper + lower) / 200;  % e.g. '100-90'→0.95
end

%% Create interpolants R0,R1,R2,C1,C2(T,SOC)
 F_R0 = scatteredInterpolant(param.T, numSOC, param.R0, 'linear', 'nearest');
        F_R1 = scatteredInterpolant(param.T, numSOC, param.R1, 'linear', 'nearest');
        F_R2 = scatteredInterpolant(param.T, numSOC, param.R2, 'linear', 'nearest');
        F_C1 = scatteredInterpolant(param.T, numSOC, param.C1, 'linear', 'nearest');
        F_C2 = scatteredInterpolant(param.T, numSOC, param.C2, 'linear', 'nearest');

%% OCV-SOC curve
%% OCV-SOC curve
SOCOCV = polyfit(SOC_OCV_AUKF.SOC, SOC_OCV_AUKF.OCV, 12);
dSOCOCV = polyder(SOCOCV);

% Create inverse OCV-SOC relationship for initial SOC estimation 9
OCV_INVSOC = polyfit(SOC_OCV_AUKF.OCV, SOC_OCV_AUKF.SOC, 8);

%% Define SOC regions for Q/R adaptation - More aggressive for better tracking
soc_regions = 1.0:-0.1:0.0;     % [1.0,0.9,...,0.0]
num_regions = numel(soc_regions) - 1;
Q_regions = cell(num_regions, 1);
R_regions = zeros(num_regions, 1);

% Much more aggressive noise parameters to track measured SOC
for i = 1:num_regions
    if i <= 1  % High SOC regions (100%-80%)
        Q_regions{i} = diag([1.3*1e-10, 2.88*1e-7, 8.96*1e-8]); R_regions(i) = 6.13*1e-3;
     elseif i <= 2  % Medium-high SOC (80%-60%)
        Q_regions{i} = diag([1.21*1e-10, 3.47*1e-11, 6.12*1e-12]); R_regions(i) = 1.08*1e-1;
     elseif i <= 3  % Medium-high SOC (80%-60%)
        Q_regions{i} = diag([4.98e-9, 3.20*1e-5, 4.43*1e-11]); R_regions(i) = 1.34*1e-4;
    elseif i <= 4  % Medium SOC (60%-40%)
        Q_regions{i} = diag([4.05*1e-11, 2.04*1e-8, 2.11*1e-10]); R_regions(i) = 1.38e-1;
    elseif i <= 5  % Medium-high SOC (80%-60%)
        Q_regions{i} = diag([2.58*1e-10, 1.24*1e-4, 4.26*1e-09]); R_regions(i) = 1.33*1e-2;
     elseif i <= 6  % Medium-high SOC (80%-60%)
        Q_regions{i} = diag([6.58e-11, 4.58*1e-8, 4.32*1e-09]); R_regions(i) = 1.22*1e-0;
     elseif i <= 7  % Medium-high SOC (80%-60%)
        Q_regions{i} = diag([3.6*1e-10, 2.2*1e-12, 1.4*1e-04]); R_regions(i) = 4.19*1e-0;
     elseif i <= 8  % Medium-high SOC (80%-60%)
        Q_regions{i} = diag([1.17e-12, 1.80*1e-7, 8.6*1e-06]); R_regions(i) = 4.47*1e-2;
    elseif i <= 9  % Medium-high SOC (80%-60%)
        Q_regions{i} = diag([8.79e-9, 4.62*1e-3, 1.26*1e-03]); R_regions(i) = 3.04*1e-1;
    else           % Very low SOC (20%-0%)
        Q_regions{i} = diag([8.42*1e-7, 3.88*1e-5, 1.62*1e-9]); R_regions(i) = 0.936;
    end
end

% Function to determine current SOC region
get_soc_region = @(soc) find(soc <= soc_regions(1:end-1) & soc > soc_regions(2:end), 1);

%% ASRUKF Parameters - More aggressive for tracking
n = 3; % State dimension [SOC, V1, V2]
alpha = 1; % Increased alpha for better tracking
beta = 2; % UKF parameter
kappa = 3 - n; % Modified kappa
lambda = alpha^2 * (n + kappa) - n; % UKF parameter

% Weights
w_m = zeros(2*n+1, 1);
w_c = zeros(2*n+1, 1);
w_m(1) = lambda/(n+lambda);
w_c(1) = w_m(1) + (1 - alpha^2 + beta);
for i = 2:2*n+1
    w_m(i) = 1/(2*(n+lambda));
    w_c(i) = w_m(i);
end

%% Initialize ASRUKF with fixed initial SOC
initial_soc_estimate = 0.8;
fprintf('Initial SOC set to: %.3f (%.1f%%)\n', initial_soc_estimate, initial_soc_estimate*100);

x_hat = [initial_soc_estimate; 0; 0]; % Initial state estimate [SOC, V1, V2]

% Initial covariance - larger uncertainty to allow adaptation
P_init = [1.25 0 0;      % Larger initial SOC uncertainty
         0 0.05 0;       % Small V1 uncertainty
         0 0 0.01];      % Small V2 uncertainty
S = chol(P_init, 'lower'); % Initial covariance square root

% Initialize with appropriate region's noise parameters
current_region = get_soc_region(initial_soc_estimate);
if isempty(current_region)
    current_region = 1; % Default to first region if outside bounds
end

SQ = chol(Q_regions{current_region}, 'lower'); % Process noise square root
SR = sqrt(R_regions(current_region)); % Measurement noise

% Store previous region for tracking changes
previous_region = current_region;

% Region change adaptation parameters
region_change_count = 0;
samples_in_region = 0;
min_samples_for_adaptation = 10; % More samples before adaptation

%% Initialize output arrays
N = numel(Current);
SOC_Estimated = zeros(N, 1);
Vt_Estimated = zeros(N, 1);
Vt_Error = zeros(N, 1);
Region_History = zeros(N, 1); % Track which region is being used

%% Main ASRUKF Loop
for k = 1:N
    T = Temperature(k); % Current temperature
    U = Current(k); % Current measurement
    SOC = x_hat(1);
    V1 = x_hat(2);
    V2 = x_hat(3);
    
    % Determine current SOC region
    current_region = get_soc_region(SOC);
    if isempty(current_region)
        % Handle edge cases
        if SOC >= 1.0
            current_region = 1;
        else
            current_region = num_regions;
        end
    end
    
    % Check if region has changed and update noise parameters
    if current_region ~= previous_region && samples_in_region >= min_samples_for_adaptation
        fprintf('SOC Region changed from %d to %d at sample %d (SOC: %.3f)\n', ...
                previous_region, current_region, k, SOC);
        
        % Update noise parameters for new region
        SQ = chol(Q_regions{current_region}, 'lower');
        SR = sqrt(R_regions(current_region));
        
        previous_region = current_region;
        region_change_count = region_change_count + 1;
        samples_in_region = 0;
    end
    
    samples_in_region = samples_in_region + 1;
    Region_History(k) = current_region;
    
    % Get battery parameters at current temperature & SOC
    R0 = F_R0(T, SOC);
    R1 = F_R1(T, SOC);
    R2 = F_R2(T, SOC);
    C1 = F_C1(T, SOC);
    C2 = F_C2(T, SOC);
    
    % Calculate time constants and discrete-time parameters
    Tau_1 = C1 * R1;
    Tau_2 = C2 * R2;
    a1 = exp(-deltaT/Tau_1);
    a2 = exp(-deltaT/Tau_2);
    b1 = R1 * (1 - exp(-deltaT/Tau_1));
    b2 = R2 * (1 - exp(-deltaT/Tau_2));
    
    % OCV using polynomial
    OCV = polyval(SOCOCV, SOC);
    
    % Calculate terminal voltage for output
    TerminalVoltage = OCV - R0*U - V1 - V2;
    
    % Calculate error
    Error_x = Vt_Actual(k) - TerminalVoltage;
    
    % Store outputs
    Vt_Estimated(k) = TerminalVoltage;
    SOC_Estimated(k) = x_hat(1);
    Vt_Error(k) = Error_x;
    
    % Skip ASRUKF update for first sample
    if k == 1
        continue;
    end
    
    % System matrices
    A = [1, 0, 0;
         0, a1, 0;
         0, 0, a2];
    
    eta_current = 1; 
    
    B = [-(eta_current * deltaT/QN);
         b1;
         b2];
    
    %% ASRUKF Steps
    
    % Step 1: Generate sigma points
    sqrt_term = sqrt(n + lambda) * S;
    chi = [x_hat, x_hat + sqrt_term, x_hat - sqrt_term];
    
    % Step 2: Time update
    chi_pred = zeros(3, 2*n+1);
    for i = 1:2*n+1
        chi_pred(:, i) = A * chi(:, i) + B * U;
        % Ensure SOC bounds
        chi_pred(1, i) = max(0.01, min(0.99, chi_pred(1, i)));
    end
    
    % Predicted state mean
    x_pred = zeros(3, 1);
    for i = 1:2*n+1
        x_pred = x_pred + w_m(i) * chi_pred(:, i);
    end
    
    % Predicted covariance square root using QR decomposition
    X_dev = chi_pred(:, 2:end) - x_pred;
    for i = 1:2*n
        X_dev(:, i) = sqrt(abs(w_c(i+1))) * X_dev(:, i);
    end
    
    [~, S_pred] = qr([X_dev, SQ]', 0);
    S_pred = S_pred';
    
    % Cholupdate for first sigma point
    if w_c(1) >= 0
        S_pred = cholupdate(S_pred, sqrt(w_c(1)) * (chi_pred(:, 1) - x_pred), '+');
    else
        S_pred = cholupdate(S_pred, sqrt(abs(w_c(1))) * (chi_pred(:, 1) - x_pred), '-');
    end
    
    % Step 3: Generate new sigma points for measurement update
    sqrt_term = sqrt(n + lambda) * S_pred;
    chi_pred = [x_pred, x_pred + sqrt_term, x_pred - sqrt_term];
    
    % Step 4: Measurement prediction
    z_pred = zeros(1, 2*n+1);
    
    for i = 1:2*n+1
        soc_i = max(0.01, min(0.99, chi_pred(1, i))); % Ensure SOC bounds
        v1_i = chi_pred(2, i);
        v2_i = chi_pred(3, i);
        
        % Get parameters for this sigma point
        R0_i = F_R0(T, soc_i);
        
        % OCV using polynomial
        OCV_i = polyval(SOCOCV, soc_i);
        
        % Terminal voltage prediction
        z_pred(i) = OCV_i - R0_i*U - v1_i - v2_i;
    end
    
    % Predicted measurement mean
    z_mean = 0;
    for i = 1:2*n+1
        z_mean = z_mean + w_m(i) * z_pred(i);
    end
    
    % Step 5: Calculate Kalman gain
    
    % Calculate measurement covariance matrix elements
    Z_dev = z_pred(2:end) - z_mean;
    for i = 1:2*n
        Z_dev(i) = sqrt(abs(w_c(i+1))) * Z_dev(i);
    end
    
    % Calculate Szz using QR decomposition  
    [~, Szz] = qr([Z_dev, SR]', 0);
    Szz = Szz(1,1); % Extract scalar result
    
    % Cholupdate for Szz with first sigma point
    if w_c(1) >= 0
        Szz = cholupdate(Szz, sqrt(abs(w_c(1))) * (z_pred(1) - z_mean), '+');
    else
        Szz = cholupdate(Szz, sqrt(abs(w_c(1))) * (z_pred(1) - z_mean), '-');
    end
    
    % Calculate cross-covariance Pxz
    Pxz = zeros(3, 1);
    for i = 1:2*n+1
        Pxz = Pxz + w_c(i) * (chi_pred(:, i) - x_pred) * (z_pred(i) - z_mean);
    end
    
    % Calculate Kalman gain
    K = Pxz / (Szz * Szz');
    
    % Step 6: State and covariance update
    innovation = Vt_Actual(k) - z_mean;
    x_hat = x_pred + K * innovation;
    
    % Ensure SOC bounds
    x_hat(1) = max(0.01, min(0.99, x_hat(1)));
    
    % Update covariance square root - Use Joseph form for stability
    try
        S = cholupdate(S_pred, K * Szz, '-');
    catch
        % If cholupdate fails, use Joseph form update
        H_jacobian = [polyval(dSOCOCV, x_hat(1)), -1, -1];
        I_KH = eye(3) - K * H_jacobian;
        P_new = I_KH * (S_pred * S_pred') * I_KH' + K * (SR^2) * K';
        % Add regularization
        P_new = P_new + 1e-10 * eye(3);
        try
            S = chol(P_new, 'lower');
        catch
            % If still fails, use SVD for positive definite matrix
            [U, Sigma, V] = svd(P_new);
            Sigma = max(Sigma, 1e-8); % Ensure positive eigenvalues
            P_new = U * Sigma * V';
            S = chol(P_new, 'lower');
        end
    end
end

%% Summary
fprintf('\nImproved Region-based ASRUKF Summary:\n');
fprintf('Initial SOC: %.1f%%\n', initial_soc_estimate*100);
fprintf('Final SOC: %.1f%%\n', x_hat(1)*100);
fprintf('Total region changes: %d\n', region_change_count);
fprintf('Final SOC region: %d (%.1f%% - %.1f%%)\n', current_region, ...
        soc_regions(current_region)*100, soc_regions(current_region+1)*100);

% Display region usage statistics
region_usage = zeros(num_regions, 1);
for i = 1:num_regions
    region_usage(i) = sum(Region_History == i);
end
fprintf('Region usage (samples per region):\n');
for i = 1:num_regions
    fprintf('  Region %d (%.0f%%-%.0f%%): %d samples\n', i, ...
            soc_regions(i)*100, soc_regions(i+1)*100, region_usage(i));
end

end