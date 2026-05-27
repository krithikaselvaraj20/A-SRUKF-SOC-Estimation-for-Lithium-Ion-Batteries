% ECM parameters (from EKF) - ORIGINAL REGION 5
R0 = 0.0012;  % 0.1 mΩ
R1 = 0.0012; % 6.21 mΩ  
R2 = 0.0005; % 2.85 mΩ
C1 = 32308.6;    % 9.9 kF
C2 = 34720.7;   % 14.6 kF
OCV = 3.288; % Initial OCV

% Time, Current, and Voltage data
t_data = [23400, 23400.5, 23401, 23401.5, 23402, 23402.5, 23403, 23403.5, ...
        23404, 23404.5, 23405, 23405.5, 23406, 23406.5, 23407, 23407.5, ...
        23408, 23408.5, 23409, 23409.5, 23410, 23410.5]; % Time (s)
    
I_data = [0, -100.0122, -100.0122, -100.0122, -100.0156, -100.0122, ...
           -100.0156, -100.0122, -100.0156, -100.0122, -100.0122, ...
           -100.0122, -100.0156, -100.0156, -100.0122, -100.0156, ...
           -100.0156, -100.0122, -100.0156, -100.0156, -100.0122, 0]; % Current (A)
      
V_data = [3.2882, 3.1589, 3.152, 3.1471, 3.1434, 3.1407, 3.1382, ...
              3.1363, 3.1344, 3.1328, 3.1313, 3.1299, 3.1287, 3.1274, ...
              3.1263, 3.1253, 3.1244, 3.1234, 3.1224, 3.1216, 3.1208, 3.2863]; % Voltage (V)
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