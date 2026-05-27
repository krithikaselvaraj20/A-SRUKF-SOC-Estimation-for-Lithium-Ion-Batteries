

% Your ECM parameters (from EKF)
R0=0.0008;
R1=0.0011;
R2=0.0005;
C1=30506.1;
C2=34084.4;
OCV =3.381; % Initial OCV

% Time, Current, and Voltage data
% Experimental Data
t_data = [300, 300.5, 301, 301.5, 302, 302.5, 303, 303.5, 304, 304.5, ...
        305, 305.5, 306, 306.5, 307, 307.5, 308, 308.5, 309, 309.5, 310, 310.5];
I_data = [0, -100.0053, -100.0088, -100.0053, -100.0088, -100.0088, -100.0088, -100.0088, ...
           -100.0088, -100.0088, -100.0088, -100.0122, -100.0088, -100.0122, -100.0088, ...
           -100.0122, -100.0088, -100.0088, -100.0088, -100.0088, -100.0122, 0];
V_data = [3.5078, 3.3807, 3.3529, 3.3306, 3.315, 3.304, 3.2959, 3.2893, 3.284, 3.2795, ...
              3.2756, 3.2722, 3.2691, 3.2664, 3.2638, 3.2614, 3.2593, 3.2572, 3.2554, 3.2535, 3.2517, 3.3253];

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