%% Read and process Excel data
data = readtable('Amptons_constant_discharge.xlsx');
Time = seconds(duration(data.Time, 'InputFormat', 'hh:mm:ss'));
Voltage = data.Voltage_V_;
Current = data.Current_A_;
Ah = data.Capacity_Ah_;
Battery_Temp_degC = data.Temperature;
SOC_measured = data.SoC___;

% Create meas struct
meas.Time = Time;
meas.Voltage = Voltage;
meas.Current = Current;
meas.Ah = Ah;
meas.Battery_Temp_degC = Battery_Temp_degC;
meas.SOC_measured = SOC_measured;

% Optional: Save this if needed for future runs
save('meas_from_excel.mat', 'meas');

% Organize the data into LiPoly structure
LiPoly.RecordingTime = meas.Time;
LiPoly.Measured_Voltage = meas.Voltage;
LiPoly.Measured_Current = meas.Current;
LiPoly.Measured_Temperature = meas.Battery_Temp_degC;
LiPoly.Measured_SOC = meas.SOC_measured;


% Reverse current sign for EKF/UKF convention
LiPoly.Measured_Current_R = -LiPoly.Measured_Current;

% Convert time to hours for plotting
LiPoly.RecordingTime_Hours = LiPoly.RecordingTime / 3600;

% Display data info
fprintf('Total samples being processed: %d\n', length(LiPoly.RecordingTime));
fprintf('Time range: %.2f to %.2f hours\n', min(LiPoly.RecordingTime_Hours), max(LiPoly.RecordingTime_Hours));
fprintf('SOC range: %.1f%% to %.1f%%\n', min(LiPoly.Measured_SOC), max(LiPoly.Measured_SOC));

% Run ASRUKF estimation
[SOC_Estimated, Vt_Estimated, Vt_Error] = i_write4_aukf_2( ...
    LiPoly.Measured_Current_R, ...
    LiPoly.Measured_Voltage, ...
    LiPoly.Measured_Temperature, ...
    LiPoly.RecordingTime);

% Store estimation results
LiPoly.Estimated_SOC = SOC_Estimated;
LiPoly.Estimated_Voltage = Vt_Estimated;
LiPoly.Voltage_Error = Vt_Error;

% Calculate performance metrics
RMSE_SOC = sqrt(mean((LiPoly.Measured_SOC - 100*SOC_Estimated).^2));
RMSE_Voltage = sqrt(mean((LiPoly.Measured_Voltage - Vt_Estimated).^2));
MAE_SOC = mean(abs(LiPoly.Measured_SOC - 100*SOC_Estimated));
MAE_Voltage = mean(abs(LiPoly.Measured_Voltage - Vt_Estimated));

fprintf('\nPerformance Metrics:\n');
fprintf('SOC RMSE: %.2f%%\n', RMSE_SOC);
fprintf('SOC MAE: %.2f%%\n', MAE_SOC);
fprintf('Voltage RMSE: %.4f V\n', RMSE_Voltage);
fprintf('Voltage MAE: %.4f V\n', MAE_Voltage);
figure;
% Plot SOC comparison
plot(LiPoly.RecordingTime_Hours, LiPoly.Measured_SOC, 'b', 'LineWidth', 1.5);
hold on;
plot(LiPoly.RecordingTime_Hours, 100*SOC_Estimated, 'r--', 'LineWidth', 1.5);
legend('Measured SOC', 'Estimated SOC ASRUKF', 'Location', 'best');
xlabel('Time [hr]');
ylabel('SOC [%]');
title('Measured vs. Estimated SOC (ASRUKF)');
grid on;


fprintf('ASRUKF SOC Estimation complete using Excel-based input.\n');

Vt_measured = LiPoly.Measured_Voltage;
Vt_estimated = LiPoly.Estimated_Voltage;

Vt_RMSE   = sqrt(mean((Vt_measured - Vt_estimated).^2));
Vt_MAE    = mean(abs(Vt_measured - Vt_estimated));
Vt_MaxAE  = max(abs(Vt_measured - Vt_estimated));

% ----- Error Metrics for SOC -----
SOC_measured = LiPoly.Measured_SOC;        % Already in percentage
SOC_estimated = 100 * LiPoly.Estimated_SOC;

SOC_RMSE   = sqrt(mean((SOC_measured - SOC_estimated).^2));
SOC_MAE    = mean(abs(SOC_measured - SOC_estimated));
SOC_MaxAE  = max(abs(SOC_measured - SOC_estimated));

% ----- Display the results -----
fprintf('\n--- Terminal Voltage Estimation Errors ---\n');
fprintf('RMSE:    %.4f V\n', Vt_RMSE);
fprintf('MAE:     %.4f V\n', Vt_MAE);
fprintf('Max AE:  %.4f V\n', Vt_MaxAE);

fprintf('\n--- SOC Estimation Errors ---\n');
fprintf('RMSE:    %.4f %%\n', SOC_RMSE);
fprintf('MAE:     %.4f %%\n', SOC_MAE);
fprintf('Max AE:  %.4f %%\n', SOC_MaxAE);
