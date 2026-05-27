% ECM parameters (from EKF) - ORIGINAL REGION 2
R0 = 0.0012;  % 0.1 mΩ
R1 = 0.0010; % 6.21 mΩ  
R2 = 0.0004; % 2.85 mΩ
C1 = 31132.8;    % 9.9 kF
C2 = 38227.7;   % 14.6 kF
OCV = 3.3253; % Initial OCV

% Time, Current, and Voltage data
t_data = [4920, 4920.5, 4921, 4921.5, 4922, 4922.5, 4923, 4923.5, 4924, 4924.5, 4925, 4925.5, 4926, 4926.5, 4927, 4927.5, 4928, 4928.5, 4929, 4929.5, 4930, 4930.5];
I_data = [0, -100.0088, -100.0088, -100.0088, -100.0122, -100.0088, -100.0122, -100.0122, -100.0122, -100.0122, -100.0088, -100.0088, -100.0122, -100.0088, -100.0088, -100.0088, -100.0088, -100.0122, -100.0088, -100.0122, -100.0088, 0];
V_data = [3.3253, 3.2062, 3.2017, 3.1983, 3.1959, 3.1938, 3.1922, 3.1905, 3.1891, 3.1878, 3.1867, 3.1855, 3.1844, 3.1834, 3.1825, 3.1815, 3.1807, 3.1799, 3.1789, 3.1781, 3.1775, 3.3261];

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