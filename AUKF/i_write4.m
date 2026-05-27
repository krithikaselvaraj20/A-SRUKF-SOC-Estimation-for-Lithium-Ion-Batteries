function [SOC_Estimated, Vt_Estimated, Vt_Error] = i_write4(Current, Vt_Actual, Temperature, RecordingTime)
% ASRUKF_SOC_Estimation - Adaptive Square Root Unscented Kalman Filter for SOC Estimation
% Adapted to use EKF parameters and real-time data
%
% Inputs:
%   Current     - Current measurements (A) [N x 1]
%   Vt_Actual   - Actual terminal voltage measurements (V) [N x 1]
%   Temperature - Temperature measurements (°C) [N x 1]
%   RecordingTime - Time vector [N x 1]
%
% Outputs:
%   SOC_Estimated - Estimated SOC values [N x 1]
%   Vt_Estimated  - Estimated terminal voltage values (V) [N x 1]
%   Vt_Error      - Terminal voltage estimation error (V) [N x 1]

%% Load Battery Parameters (same as EKF)
load 'BatteryModel.mat'; % Load the battery parameters 
load 'SOC-OCV.mat'; % Load the SOC-OCV curve

%% Battery Parameters (from EKF code)
SOC_Init = 0.7; % Initial SOC
QN = 4.81 * 3600; % Nominal capacity (Ah to Amp-seconds)
eta = 1.0; % Coulombic efficiency

if length(RecordingTime) > 1
    time_diffs = RecordingTime(2:end) - RecordingTime(1:end-1);
    deltaT = median(time_diffs); % Use median for robustness
    fprintf('Calculated sampling time: %.4f seconds\n', deltaT);
else
    deltaT = 9.85;
    fprintf('Using default sampling time: %.4f seconds\n', deltaT);
end

% Initialize scatteredInterpolant functions for battery parameters
F_R0 = scatteredInterpolant(param.T, param.SOC, param.R0);
F_R1 = scatteredInterpolant(param.T, param.SOC, param.R1);
F_R2 = scatteredInterpolant(param.T, param.SOC, param.R2);
F_C1 = scatteredInterpolant(param.T, param.SOC, param.C1);
F_C2 = scatteredInterpolant(param.T, param.SOC, param.C2);

% SOC-OCV relationship using polynomial (from EKF code)
SOCOCV = polyfit(SOC_OCV.SOC, SOC_OCV.OCV, 6); % 8th order polynomial
dSOCOCV = polyder(SOCOCV); % derivative for Jacobian

%% ASRUKF Parameters
n = 3; % State dimension [SOC, V1, V2]
alpha = 0.01; % Reduced alpha for better performance
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

%% Initialize ASRUKF
% Initialize SOC using first voltage measurement for better starting point
initial_voltage = Vt_Actual(1);
initial_current = Current(1);
initial_temp = Temperature(1);

% Estimate initial SOC from first measurement
% Use a simple voltage-based initialization with wider search range

initial_soc_estimate = 0.8;

fprintf('Initial SOC estimate from voltage: %.3f\n', initial_soc_estimate);

x_hat = [initial_soc_estimate; 0; 0]; % Initial state estimate [SOC, V1, V2]

% Adjusted initial covariance - increased SOC uncertainty to allow better adaptation
P_init = [0.25 0 0;      % Higher initial SOC uncertainty
         0 0.05 0;
         0 0 0.05];
S = chol(P_init, 'lower'); % Initial covariance square root

% Adjusted process noise - higher for SOC to allow better tracking
Q_init = [5e-4 0 0;      % Higher SOC process noise
         0 1.0e-6 0;
         0 0 1.0e-6];
SQ = chol(Q_init, 'lower'); % Process noise square root
SR = sqrt(1.0e-4); % Measurement noise (reduced from original)

% Minimum noise bounds
SQ_min = diag([1.0e-6, 1.0e-8, 1.0e-8]);
SR_min = 1.0e-5;


% Adaptive parameters
b = 0.95; % Increased forgetting factor for smoother adaptation
M = 10; % Increased window size
ev_history = zeros(M, 1);
innovation_count = 0;

%% Initialize output arrays
ik = length(Current);
SOC_Estimated = zeros(ik, 1);
Vt_Estimated = zeros(ik, 1);
Vt_Error = zeros(ik, 1);

%% Main ASRUKF Loop
for k = 1:ik
    T = Temperature(k); % Current temperature
    U = Current(k); % Current measurement
    SOC = x_hat(1);
    V1 = x_hat(2);
    V2 = x_hat(3);
    
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
        P_new = P_new + 1e-8 * eye(3);
        S = chol(P_new, 'lower');
    end
    
    % Step 7: Adaptive noise estimation
    innovation_count = innovation_count + 1;
    
    % Update residual history
    if innovation_count <= M
        ev_history(innovation_count) = innovation^2;
    else
        ev_history = [ev_history(2:end); innovation^2];
    end
    
    % Only adapt after sufficient samples
    if innovation_count >= M
        F_k = mean(ev_history);
        
        % Conservative adaptation
           d_k = (1 - b^k)/(1 - b^(k+1));  % Reduced adaptation rate for stability
        
        % Only adapt if innovation is reasonable
        %if F_k > 1e-8 && F_k < 1e-2
            % Process noise adaptation - more aggressive for SOC
            Q_factor = sqrt(F_k);
            SQ_new = diag([Q_factor*1e-3 Q_factor*1e-5, Q_factor*1e-5]) + SQ_min;
            SQ = (1 - d_k) * SQ + d_k * SQ_new;
            
            % Measurement noise adaptation  
            SR_new = sqrt(F_k)+ SR_min;
            SR = (1 - d_k) * SR + d_k * SR_new;
        %end
    end
end

end