clc; clear; close all;
% Load the battery data
load('06-03-19_09.46 825_LA92_0degC_Turnigy_Graphene.mat');
% Organize the data into LiPoly structure
LiPoly.RecordingTime = meas.Time;
LiPoly.Measured_Voltage = meas.Voltage;
LiPoly.Measured_Current = meas.Current;
LiPoly.Measured_Temperature = meas.Battery_Temp_degC;
% Battery nominal capacity in Ah
nominalCap = 4.81; % Battery capacity in Ah taken from data
% Calculate the SOC using Coulomb Counting for comparison
LiPoly.Measured_SOC = (nominalCap + meas.Ah) .* 100 ./ nominalCap;
% Resample input data (downsample by factor of 10)
LiPoly.RecordingTime = LiPoly.RecordingTime(1:10:end);
LiPoly.Measured_Voltage = LiPoly.Measured_Voltage(1:10:end);
LiPoly.Measured_Current = LiPoly.Measured_Current(1:10:end);
LiPoly.Measured_Temperature = LiPoly.Measured_Temperature(1:10:end);
LiPoly.Measured_SOC = LiPoly.Measured_SOC(1:10:end);
% Current Definition: (+) Discharging, (-) Charging
% Reverse the current sign to match EKF convention
LiPoly.Measured_Current_R = - LiPoly.Measured_Current;
% Converting seconds to hours
LiPoly.RecordingTime_Hours = LiPoly.RecordingTime/3600;
% Call ASRUKF SOC Estimation function with same parameters as EKF
% Make sure your ASRUKF function signature matches: (Current, Voltage, Temperature)
[SOC_Estimated, Vt_Estimated, Vt_Error] = i_write4(LiPoly.Measured_Current_R, ...
 LiPoly.Measured_Voltage, ...
 LiPoly.Measured_Temperature,...
 LiPoly.RecordingTime);
% Store the results
LiPoly.Estimated_SOC = SOC_Estimated;
LiPoly.Estimated_Voltage = Vt_Estimated;
LiPoly.Voltage_Error = Vt_Error;

% ----- Error Metrics for Terminal Voltage -----
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

% Create the same plot as your EKF version
figure;
plot(LiPoly.RecordingTime_Hours, LiPoly.Measured_Voltage);
hold on;
plot(LiPoly.RecordingTime_Hours, Vt_Estimated);
hold off;
legend('Measured', 'Estimated ASRUKF');
ylabel('Terminal Voltage[V]');
xlabel('Time[hr]');
title('Measured vs. Estimated Terminal Voltage (V) at 0 Deg C');
grid minor;

% Add SOC comparison plot
figure;
plot(LiPoly.RecordingTime_Hours, LiPoly.Measured_SOC);
hold on;
plot(LiPoly.RecordingTime_Hours, 100*SOC_Estimated);
hold off;
legend('Measured SOC', 'Estimated SOC ASRUKF');
ylabel('State of Charge [%]');
xlabel('Time[hr]');
title('Measured vs. Estimated SOC (%) at 0 Deg C');
grid minor;

% Display results
fprintf('ASRUKF SOC Estimation completed.\n');
fprintf('Data points processed: %d\n', length(LiPoly.RecordingTime));