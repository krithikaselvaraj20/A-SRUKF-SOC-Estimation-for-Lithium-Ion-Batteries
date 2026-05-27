%% EKF-based RC Parameter Estimation for HPPC Data (with region-specific OCV)
rng(42);  % Sets the random number seed, consistent every time

clear; clc; close all;
%% Configuration
num_regions = 9;
results = cell(num_regions, 1);
%% Process each region
for region = 1:num_regions
    fprintf('Processing Region %d...\n', region);
    filename = sprintf('discharge_pts_region%d.xlsx', region);
    try
        data = readtable(filename);
        if width(data) >= 3
            time = table2array(data(:,1));
            current = table2array(data(:,2));
            voltage = table2array(data(:,3));
        else
            error('Excel file must have at least 3 columns: Time, Current, Voltage');
        end

        raw = readmatrix(filename);
        ocv_value = raw(1,3);

        [estimated_params, estimated_states, innovations] = ...
            ekf_rc_estimation_modified(time, voltage, current, region, ocv_value);
        results{region} = struct(...
            'region', region, 'R0', estimated_params.R0, 'R1', estimated_params.R1, ...
            'C1', estimated_params.C1, 'R2', estimated_params.R2, 'C2', estimated_params.C2, ...
            'tau1', estimated_params.R1 * estimated_params.C1, ...
            'tau2', estimated_params.R2 * estimated_params.C2, ...
            'rmse', estimated_params.rmse, 'states', estimated_states, ...
            'innovations', innovations);
        fprintf('Region %d completed successfully\n', region);
    catch ME
        fprintf('Error processing Region %d: %s\n', region, ME.message);
        results{region} = struct('region', region, 'error', ME.message);
    end
end
%% Display Results Summary
fprintf('\n=== RC Parameter Estimation Results ===\n');
fprintf('Region\tR0(Ω)\t\tR1(Ω)\t\tC1(F)\t\tR2(Ω)\t\tC2(F)\t\tτ1(s)\t\tτ2(s)\n');

fprintf('------\t------\t\t------\t\t------\t\t------\t\t------\t\t------\t\t------\t\t\n');
for region = 1:num_regions
    if isfield(results{region}, 'R0')
        r = results{region};
        fprintf('%d\t%.4f\t\t%.4f\t\t%.1f\t\t%.4f\t\t%.1f\t\t%.2f\t\t%.2f\n', ...
            r.region, r.R0, r.R1, r.C1, r.R2, r.C2, r.tau1, r.tau2);

    else
        fprintf('%d\tError in processing\n', region);
    end
end
%% Save results
save('hppc_rc_parameters_modified.mat', 'results');
fprintf('\nResults saved to hppc_rc_parameters_modified.mat\n');

%% EKF RC Estimation Function
function [estimated_params, estimated_states, innovations] = ekf_rc_estimation_modified(time, voltage, current, region, OCV)
    dt = mean(diff(time)); N = length(time);
    if region == 1
        R0_init = 0.001 * (1 + 0.01*randn());  
    else
        R0_init = 0.02 * (1 + 0.1*randn());
    end
    if region == 2 || region == 3
        R0_init = 0.0009 * (1 + 0.02*randn());  
    else
        R0_init = 0.02 * (1 + 0.1*randn());
    end
    R1_init = 0.001 * (1 + 0.1*randn());
    C1_init = 30000 * (1 + 0.05*randn());
    R2_init = 0.0005 * (1 + 0.1*randn());
    C2_init = 35000 * (1 + 0.05*randn());
    
    x = [0; 0; R0_init; R1_init; C1_init; R2_init; C2_init];
    P = diag([0.01, 0.01, 1e-6, 1e-8, 1e4, 1e-8, 1e4]);
    Q = diag([1e-4, 1e-4, 1e-6, 1e-9, 5e3, 1e-9, 5e3]) * (1 + 0.05*region);
    Q(5,5) = 1e5; 
    Q(7,7) = 1e5;  
    if region == 1
        Q = diag([1e-4, 1e-4, 1e-7, 1e-9, 3e3, 1e-9, 3e3]); % Larger Q(3,3)
    end
    if region == 2 || region == 3
        Q = diag([1e-4, 1e-4, 1e-7, 1e-9, 4e3, 1e-9, 4e3]);  % Conservative Q for stable C1/C2  
    end
    R_noise = get_measurement_noise_variance_modified(region);

    estimated_states = zeros(7, N);
    innovations = zeros(1, N);
    voltage_est = zeros(N, 1);

    for k = 1:N
        % 1.state prediction
        [x_pred, F] = predict_state_modified(x, current(k), dt);

        % 2.covariance prediction
        P_pred = F * P * F' + Q;
        
        % 3.measurement
        [z_pred, H] = predict_measurement_modified(x_pred, current(k), OCV);

        % 4.innovation
        innovation = voltage(k) - z_pred; % Difference between actual and predicted
        
        % 5.kalmann gain
        S = H * P_pred * H' + R_noise;
        K = P_pred * H' / S;

        % 6.updation
        x = x_pred + K * innovation;
        P = (eye(7) - K * H) * P_pred;


        % Bound parameters
        if region == 1
            x(3) = max(min(x(3), 0.1), 0.0008); 
        elseif region == 7 || region == 8 || region == 9
            x(3) = max(min(x(3), 0.1), 0.00145); 
        else
            x(3) = max(min(x(3), 0.1), 0.0012); 
        end
        x(4) = max(min(x(4), 0.01), 1e-6);       % R1
        x(5) = max(min(x(5), 50000), 20000);     % C1
        x(6) = max(min(x(6), 0.005), 1e-6);      % R2
        x(7) = max(min(x(7), 50000), 20000);     % C2
        estimated_states(:, k) = x;
        innovations(k) = innovation;
        voltage_est(k) = z_pred + innovation;
    end
    estimated_params = struct('R0', x(3), 'R1', x(4), 'C1', x(5), ...
                              'R2', x(6), 'C2', x(7), ...
                              'rmse', sqrt(mean((voltage - voltage_est).^2)));
end
%% State Transition Function
function [x_pred, F] = predict_state_modified(x, current, dt)
    V_RC1 = x(1); V_RC2 = x(2);
    R1 = x(4); C1 = x(5); R2 = x(6); C2 = x(7);
    tau1 = R1 * C1; tau2 = R2 * C2;
    exp1 = exp(-dt / tau1); exp2 = exp(-dt / tau2);
    V_RC1_pred = exp1 * V_RC1 + R1 * (1 - exp1) * current;
    V_RC2_pred = exp2 * V_RC2 + R2 * (1 - exp2) * current;
    x_pred = [V_RC1_pred; V_RC2_pred; x(3:end)];
    F = eye(7);
    F(1,1) = exp1; F(1,4) = (1 - exp1) * current;
    F(2,2) = exp2; F(2,6) = (1 - exp2) * current;
end

%% Measurement Function
function [z_pred, H] = predict_measurement_modified(x, current, OCV)
    V_RC1 = x(1); V_RC2 = x(2); R0 = x(3);
    z_pred = OCV - V_RC1 - V_RC2 - R0 * current;
    H = [-1, -1, -current, 0, 0, 0, 0];
end
%% Measurement Noise Function
function R_noise = get_measurement_noise_variance_modified(region)
    R_base = 1e-4;
    scaling_factors = [1.0, 1.4, 0.9, 1.3, 1.2, 0.8, 1.4, 1.5, 1.6];
    R_noise = R_base * scaling_factors(region);
end  
