

% Your ECM parameters (from EKF)
R0=0.0014;
R1=0.0011;
R2=0.0005;
C1=29204.8;
C2=35271.3;
OCV =3.284; % Initial OCV

% Time, Current, and Voltage data
t_data = [32640, 32640.5, 32641, 32641.5, 32642, 32642.5, 32643, 32643.5, 32644, 32644.5, ...
          32645, 32645.5, 32646, 32646.5, 32647, 32647.5, 32648, 32648.5, 32649, 32649.5, ...
          32650, 32650.5];
I_data = [0, -100.0053, -100.0088, -100.0053, -100.0053, -100.0053, -100.0053, -100.0019, -100.0053, -100.0053, -100.0053, -100.0053, -100.0053, -100.0053, -100.0019, -100.0053, -100.0053, -100.0053, -100.0019, -100.0053, -100.0053, 0];

V_data = [3.2845, 3.1499, 3.1410, 3.1347, 3.1302, 3.1265, 3.1234, 3.1208, 3.1187, 3.1166, ...
          3.1148, 3.1131, 3.1116, 3.1102, 3.1089, 3.1076, 3.1064, 3.1053, 3.1042, 3.1032, ...
          3.1023, 3.2695];
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