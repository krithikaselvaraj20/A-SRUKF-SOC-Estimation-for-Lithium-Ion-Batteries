

% Your ECM parameters (from EKF)
R0=0.0014;
R1=0.0011;
R2=0.0006;
C1=31794.4;
C2=35678.5;
OCV =3.269; % Initial OCV

% Time, Current, and Voltage data
t_data = [37260, 37260.5, 37261, 37261.5, 37262, 37262.5, 37263, 37263.5, ...
        37264, 37264.5, 37265, 37265.5, 37266, 37266.5, 37267, 37267.5, ...
        37268, 37268.5, 37269, 37269.5, 37270]; % Time (s)
    
I_data = [0, -100.0019, -100.0019, -100.0053, -100.0053, -100.0053, ...
           -100.0019, -100.0019, -100.0019, -100.0019, -100.0053, ...
           -100.0019, -100.0019, -100.0019, -100.0053, -100.0019, ...
           -100.0053, -100.0019, -100.0019, -100.0053, -100.0019]; % Current (A)
      
V_data = [3.2695, 3.1347, 3.1248, 3.1177, 3.1126, 3.1084, 3.1048, ...
              3.1019, 3.0993, 3.0969, 3.0948, 3.0929, 3.0911, 3.0895, ...
              3.088, 3.0866, 3.0853, 3.0842, 3.0829, 3.0818, 3.0808]; % Voltage (V)
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