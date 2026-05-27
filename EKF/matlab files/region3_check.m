% ECM parameters (from EKF) - ORIGINAL REGION 2
R0 = 0.0012;  % 0.1 mΩ
R1 = 0.0009; % 6.21 mΩ  
R2 = 0.0005; % 2.85 mΩ
C1 = 28978.6;    % 9.9 kF
C2 = 34895.8;   % 14.6 kF
OCV = 3.326; % Initial OCV

% Time, Current, and Voltage data
t_data = [9540, 9540.5, 9541, 9541.5, 9542, 9542.5, 9543, 9543.5, 9544, 9544.5, ...
        9545, 9545.5, 9546, 9546.5, 9547, 9547.5, 9548, 9548.5, 9549, 9549.5, 9550, 9550.5];
    
I_data = [0, -100.0053, -100.0053, -100.0088, -100.0088, -100.0088, -100.0088, -100.0088, ...
           -100.0088, -100.0088, -100.0088, -100.0088, -100.0088, -100.0088, -100.0088, ...
           -100.0122, -100.0122, -100.0088, -100.0088, -100.0088, -100.0088, 0];
      
V_data = [3.3261, 3.2051, 3.1997, 3.196, 3.1931, 3.1909, 3.1889, 3.1872, 3.1855, ...
              3.1841, 3.1828, 3.1815, 3.1804, 3.1792, 3.1783, 3.1771, 3.1762, 3.1752, ...
              3.1742, 3.1734, 3.1726, 3.3263];

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