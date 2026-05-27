% ECM parameters (from EKF) - ORIGINAL REGION 5
R0 = 0.0014;  % 0.1 mΩ
R1 = 0.0010 % 6.21 mΩ  
R2 = 0.0005; % 2.85 mΩ
C1 = 28028;    % 9.9 kF
C2 = 33459.9;   % 14.6 kF
OCV = 3.286; % Initial OCV

% Time, Current, and Voltage data
t_data = [28020, 28020.5, 28021, 28021.5, 28022, 28022.5, 28023, 28023.5, ...
        28024, 28024.5, 28025, 28025.5, 28026, 28026.5, 28027, 28027.5, ...
        28028, 28028.5, 28029, 28029.5, 28030, 28030.5]; % Time (s)
    
I_data = [0, -100.0122, -100.0122, -100.0122, -100.0122, -100.0156, ...
           -100.0156, -100.0122, -100.0122, -100.0122, -100.0122, ...
           -100.0122, -100.0122, -100.0122, -100.0156, -100.0122, ...
           -100.0156, -100.0122, -100.0122, -100.0122, -100.0122, 0]; % Current (A)
      
V_data = [3.2863, 3.1521, 3.1445, 3.1389, 3.1347, 3.1315, 3.1289, ...
              3.1266, 3.1245, 3.1228, 3.121, 3.1195, 3.1181, 3.1169, ...
              3.1156, 3.1145, 3.1134, 3.1124, 3.1115, 3.1105, 3.1097, 3.2845]; % Voltage (V)
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