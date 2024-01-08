% Task 1: Data preprocessing and visualisation

% Load the dataset
load('Slow_Contraction.mat');
ConversionFactor = 0.02;
Gravity = 9.81; % Acceleration due to gravity g

% 1.1 Convertion of force signal to Newtons 
Force_N_Slow = ref_signal * ConversionFactor * Gravity;
time_vector_slow = (0:1/fsamp:(length(ref_signal)-1)/fsamp);

% 1.2 Plot of force signal in Newtons
figure;
plot(time_vector_slow, Force_N_Slow);
xlabel('Time (s)');
ylabel('Force (N)');
title('Force Signal - Slow Contractions');

% 1.3 Filtering the force signal
cutoff_frequency_slow = 10; % Set cutoff frequency in Hz
[b_slow, a_slow] = butter(4, cutoff_frequency_slow / (fsamp / 2), 'low');
Filtered_Force_N_Slow = filtfilt(b_slow, a_slow, Force_N_Slow);

% Plot of unfiltered and filtered signals
figure;
plot(time_vector_slow, Force_N_Slow, 'b');
hold on;
plot(time_vector_slow, Filtered_Force_N_Slow, 'r');
xlabel('Time (s)');
ylabel('Force (N)');
title('Force Signal with and without Filtering - Slow Contractions');
legend('Unfiltered', 'Filtered');

% 1.4 Plot of exemplary channel of the EMG data with force signal
channel_slow = 1;
figure;
yyaxis left;
plot(time_vector_slow, SIG{channel_slow});
ylabel('EMG Signal (mV)');

% Plot of selected channel of the EMG data and the force signal
yyaxis right;
plot(time_vector_slow, Force_N_Slow, 'r--');
ylabel('Force (N)');
xlabel('Time (s)');
title(['EMG Signal and Force Signal for Channel ' num2str(channel_slow)]);

% Task 2. Force steadiness
% 2.1 Force Steadiness and coefficient of variation (CV)
plateau_start_time_slow = 10; % seconds
plateau_end_time_slow = 20; % seconds
force_plateau_slow = Force_N_Slow(time_vector_slow >= plateau_start_time_slow & time_vector_slow <= plateau_end_time_slow);
% Calculating CV for Slow Contractions
CV_slow = std(force_plateau_slow) / mean(force_plateau_slow) * 100;
fprintf('Coefficient of Variation (CV) during plateau phase for Slow Contractions: %.2f%%\n', CV_slow);

% Task 3: Rate of Force Development for Rapid Contractions

% Task 3.1 : 
% Load the dataset 
load('Rapid_Contractions.mat');
ConversionFactor_rapid = 0.02;
Gravity_rapid = 9.81; % Acceleration due to gravity g

% Converting force signal to Newtons
Force_N_Rapid = ref_signal * ConversionFactor_rapid * Gravity_rapid;

% Obtain time vector
time_vector_rapid = (0:1/fsamp:(length(ref_signal)-1)/fsamp);

% Plotting force signal for Rapid Contractions
figure;
plot(time_vector_rapid, Force_N_Rapid);
xlabel('Time (s)');
ylabel('Force (N)');
title('Force Signal - Rapid Contractions');

% Task 3.2 :

% Onset times for each contraction
total_contractions_rapid = 5;
onset_times_rapid = linspace(5, 25, total_contractions_rapid);

rfd_window_rapid = 0.5; % 5 seconds
rfd_values_rapid = zeros(1, length(onset_times_rapid));

for i = 1:length(onset_times_rapid)
    window_indices_rapid = time_vector_rapid >= onset_times_rapid(i) & time_vector_rapid <= (onset_times_rapid(i) + rfd_window_rapid);
    rfd_values_rapid(i) = mean(diff(Force_N_Rapid(window_indices_rapid))) / mean(diff(time_vector_rapid(window_indices_rapid)));
end

% Visualizing RFD values for Rapid Contractions
figure;
bar(rfd_values_rapid);
xlabel('Contractions');
ylabel('Rate of Force Development (N/s)');
title('Rate of Force Development for Each Contraction - Rapid Contractions');


% Task 4.1: Convert SIG cell array to a 2D array and compute average across channels

% Function to normalize data
normalize_data = @(data) (data - mean(data)) / std(data);

% Slow Contractions
data_size_slow = size(SIG{1});
emg_data_slow = zeros(length(SIG), data_size_slow(2));

% Convert SIG cell array to a 2D array
for i = 1:length(SIG)
    current_size = size(SIG{i});
    emg_data_slow(i, :) = [SIG{i}, zeros(1, data_size_slow(2) - current_size(2))];
end

% Compute average across channels
average_emg_slow = mean(emg_data_slow, 1);

% Rapid Contractions
data_size_rapid = size(SIG{1});
emg_data_rapid = zeros(length(SIG), data_size_rapid(2));

% Assume the same size for rapid contractions
for i = 1:length(SIG)
    emg_data_rapid(i, :) = [SIG{i}(1:current_size(2)), zeros(1, data_size_rapid(2) - current_size(2))];

end

% Compute average across channels
average_emg_rapid = mean(emg_data_rapid, 1);



% Task 4.2: Compute the RMS of the average as a moving average with a window length of 200 ms

% RMS window length in samples
rms_window_length = round(0.2 * fsamp);

% Compute RMS for both datasets
rms_slow = rms(movmean(average_emg_slow .^ 2, rms_window_length));
rms_rapid = rms(movmean(average_emg_rapid .^ 2, rms_window_length));

% Task 4.3: Visualize the correlation between RMS and the respective force signal


% Slow Contractions
normalized_rms_slow = normalize_data(rms_slow);
normalized_force_slow = normalize_data(Force_N_Slow(1:length(rms_slow))); % Ensure lengths match

% Rapid Contractions
normalized_rms_rapid = normalize_data(rms_rapid);
normalized_force_rapid = normalize_data(Force_N_Rapid(1:length(rms_rapid))); % Ensure lengths match


% Scatter plot
figure;
subplot(2, 1, 1);
scatter(normalized_rms_slow, normalized_force_slow);
xlabel('Normalized RMS - Slow Contractions');
ylabel('Normalized Force');
title('Correlation between RMS and Force - Slow Contractions');

subplot(2, 1, 2);
scatter(normalized_rms_rapid, normalized_force_rapid);
xlabel('Normalized RMS - Rapid Contractions');
ylabel('Normalized Force');
title('Correlation between RMS and Force - Rapid Contractions');

% Compute and report the correlation coefficient R
correlation_slow = corrcoef(normalized_rms_slow, normalized_force_slow, 'Rows', 'complete');
disp('Size of correlation_slow:');
disp(size(correlation_slow));

% Check if correlation_slow is a scalar
if isscalar(correlation_slow)
    correlation_slow = correlation_slow(1);
end

correlation_rapid = corrcoef(normalized_rms_rapid, normalized_force_rapid, 'Rows', 'complete');
disp('Size of correlation_rapid:');
disp(size(correlation_rapid));

% Check if correlation_rapid is a scalar
if isscalar(correlation_rapid)
    correlation_rapid = correlation_rapid(1);
end

fprintf('Correlation coefficient R - Slow Contractions: %.2f\n', correlation_slow);
fprintf('Correlation coefficient R - Rapid Contractions: %.2f\n', correlation_rapid);



