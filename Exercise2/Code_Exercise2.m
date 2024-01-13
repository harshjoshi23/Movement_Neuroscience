

% Task 1: Spike Trains
% Add the MATLAB function directory to the search path
addpath('./Data&Codes/');

% Load the dataset for analysis
load('./Data&Codes/Slow_Contraction.mat');

% Establish a time vector for plotting if it isn't already defined
if ~exist('time_vector_slow', 'var')
    time_vector_slow = (0:length(ref_signal)-1) / fsamp;
end

% Initialize the matrix to log motor unit (MU) spike train occurrences
numMUs = length(MUPulses);  % Total number of Motor Units to consider
signalLength = length(ref_signal);  % Duration of the EMG signal in samples
firingMatrix = false(numMUs, signalLength);  % Begin with a matrix of zeros

% Fill the matrix with ones at instances where a spike is fired by an MU
for i = 1:numMUs
    firingMatrix(i, MUPulses{i}) = true;  % Mark spikes with a logical 'true'
end

% Visualizing the MU spike trains
figure;
subplot(2, 1, 1);  % Allocate a subplot for spike trains
plotSpikeRaster(logical(firingMatrix), 'PlotType', 'vertline', 'VertSpikeHeight', 0.9);
title('Motor Unit Spike Trains');

% Presenting the corresponding force signal
subplot(2, 1, 2);  % Allocate a subplot for the force signal
plot(time_vector_slow, ref_signal);  % Plot the force signal assumed to be in Newtons
xlabel('Time (s)');
ylabel('Force (N)');
title('Reference Force Signal');

% Task 2: Spike Triggered Averaging (STA)
% Determining the STA window
STA_window = 0.100; % A window of 100 milliseconds for STA computation

% Perform STA to correlate the spike trains with the continuous EMG signal
STA_result = spikeTriggeredAveraging(SIG, MUPulses, STA_window, fsamp);

% Selecting a specific motor unit for detailed examination
selectedMU = 1; % Example Motor Unit index for demonstration

% Define the electrode grid layout
nRows = 13;
nCols = 5;
totalPlots = nRows * nCols - 1;  % Consider the grid with a missing electrode

% Plot the STA results for the selected motor unit across all channels
figure(2);  % Initiate a new figure for displaying STA results
for i = 1:totalPlots
    if i == (nRows * nCols)  % Bypass the last subplot for the absent electrode
        continue
    end
    subplot(nRows, nCols, i + (i > (nRows * (nCols - 1))));  % Fit plots to the grid
    channelData = STA_result{selectedMU}{i};
    if ~isempty(channelData)
        plot(channelData);
        title(['Channel ' num2str(i)]);
    else
        plot(0);  % Insert a placeholder for absent or empty data
    end
end

MUAP_p2p_values = cell(size(STA_result));
for i = 1:numel(STA_result)
    mu = STA_result{i};
    for j = 1:numel(mu)
        ch = mu{j};
        if isempty(ch) || (numel(ch) == 1 && ch == 0)
            MUAP_p2p_values{i}{j} = NaN;
        else
            MUAP_p2p_values{i}{j} = peak2peak(ch);
        end
    end
end


% Converting the nested cell array to a numeric matrix
MUAP_p2p_matrix = cell2mat(cellfun(@(c) cell2mat(c), MUAP_p2p_values, 'UniformOutput', false));

% Debugging: Display the dimensions of the peak-to-peak matrix
disp('Size of MUAP_p2p_matrix:');
disp(size(MUAP_p2p_matrix));

% Visualize the peak-to-peak amplitudes to observe MU recruitment patterns
figure;
plot(MUAP_p2p_matrix);
xlabel('Motor Unit Index');
ylabel('Peak-to-Peak Amplitude');
title('Peak-to-Peak Amplitude of MUAPs');

% Task 4: Motor Unit Locations
% Analyzing the spatial location of motor units via MUAP shapes
selectedMU = 1; % Selecting a motor unit for location analysis

% Preparing the grid to represent the electrode layout
RMS_grid = NaN(nRows, nCols);

% Calculating the root mean square (RMS) for each channel of the selected motor unit
MUAP_RMS = zeros(1, numel(STA_result{selectedMU}));
for i = 1:numel(STA_result{selectedMU})
    channelData = STA_result{selectedMU}{i};
    if isempty(channelData)
        MUAP_RMS(i) = NaN;
    else
        MUAP_RMS(i) = rms(channelData);
    end
end

% Assigning RMS values to the grid
index = 1; % Starting index for RMS values
for row = 1:nRows
    for col = 1:nCols
        if ~(row == nRows && col == 1) % Avoiding the missing electrode spot
            RMS_grid(row, col) = MUAP_RMS(index);
            index = index + 1;
        end
    end
end

% Creating a heatmap to illustrate the MUAP RMS distribution
figure;
imagesc(RMS_grid);
colorbar; % Including a color bar for scale reference
title(sprintf('Heatmap of MUAP RMS values for Motor Unit %d', selectedMU));



