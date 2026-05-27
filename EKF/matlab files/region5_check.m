% ECM parameters (from EKF) - ORIGINAL REGION 5
R0 = 0.0012;  % 0.1 mΩ
R1 = 0.0010; % 6.21 mΩ  
R2 = 0.0004; % 2.85 mΩ
C1 = 30410.7;    % 9.9 kF
C2 = 35707;   % 14.6 kF
OCV = 3.3; % Initial OCV

% Time, Current, and Voltage data
t_data = [18780, 18780.5, 18781, 18781.5, 18782, 18782.5, 18783, 18783.5, ...
        18784, 18784.5, 18785, 18785.5, 18786, 18786.5, 18787, 18787.5, ...
        18788, 18788.5, 18789, 18789.5, 18790, 18790.5]; % Time (s)
    
I_data = [0, -100.0122, -100.0122, -100.0122, -100.0122, -100.0122, ...
           -100.0122, -100.0122, -100.0122, -100.0156, -100.0122, ...
           -100.0122, -100.0122, -100.0122, -100.0156, -100.0122, ...
           -100.0122, -100.0122, -100.0122, -100.0156, -100.0122, 0]; % Current (A)
      
V_data = [3.2998, 3.1741, 3.1678, 3.1633, 3.1599, 3.1573, 3.155, ...
              3.1531, 3.1513, 3.1497, 3.1483, 3.1468, 3.1455, 3.1444, ...
              3.1432, 3.1421, 3.1412, 3.1402, 3.1392, 3.1384, 3.1374, 3.2882]; % Voltage (V)
% Initialize simulation
n_points = length(t_data);
t_sim = t_data - t_data(1); % Normalize time to start at 0

% Initialize state variables for RC circuits
V1 = 0; % Voltage across C1
V2 = 0; % Voltage across C2
V_sim = zeros(1, n_points);

% Simulate ECM response with time-varying current
for i = 1:n_points
    if i == 1
        % Initial condition
        V_sim(i) = OCV;
        V1 = 0;
        V2 = 0;
    else
        % Time step
        dt_actual = t_data(i) - t_data(i-1);
        I_current = I_data(i);
        
        % Update RC circuit voltages using exponential decay
        % CORRECTED: Remove abs() from RC equations
        exp_factor1 = exp(-dt_actual/(R1*C1));
        exp_factor2 = exp(-dt_actual/(R2*C2));
        
        V1 = V1 * exp_factor1 - I_current * R1 * (1 - exp_factor1);
        V2 = V2 * exp_factor2 - I_current * R2 * (1 - exp_factor2);
        
        % Calculate terminal voltage
        % CORRECTED: Remove abs() from R0 term
        V_sim(i) = OCV + I_current * R0 - V1 - V2;
    end
end

% Plot comparison
figure;
plot(t_data, V_data, 'b-', 'LineWidth', 2, 'DisplayName', 'HPPC Data');
hold on;
plot(t_data, V_sim, 'r--', 'LineWidth', 2, 'DisplayName', 'ECM Simulation (Corrected)');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend;
title('ECM vs HPPC Data Validation - CORRECTED');
grid on;

% Calculate and display error metrics
rmse = sqrt(mean((V_data - V_sim).^2));
mae = mean(abs(V_data - V_sim));
max_error = max(abs(V_data - V_sim));

fprintf('Corrected Model Validation Metrics:\n');
fprintf('RMSE: %.6f V (%.2f mV)\n', rmse, rmse*1000);
fprintf('MAE: %.6f V (%.2f mV)\n', mae, mae*1000);
fprintf('Max Error: %.6f V (%.2f mV)\n', max_error, max_error*1000);