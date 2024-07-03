close all;   % Close all open figure windows
clear;       % Clear workspace

% Algorithm Parameters
movingAverageWindowSpan = 13;             % Moving average window span for signal denoising
stringencyParameter = 15;                 % Parameter for burst events detection
timeVecTotal = linspace(0, 60, 360);      % Time steps vector

% Simulation parameters
simType = 3;                                     % Simulation type 
numTracks = 1000;                                % Number of tracks to simulate
t = linspace(0, 60, 360);                        % Time vector for simulation
noiseStd = 40;                                   % Standard deviation of noise
expConst = 0.0005;                               % Exponential constant for photobleaching

% generateSim function generates simulated trajectories
% Function recieces as input: 
% numTracks - number of tracks to generate
% trackLength - the number of points in a track
% simType - 1 for constant, 2 for slope, 3 for bursts 
% expConst - constant for photobleaching exponential component
% noiseStd - standard deviation of the gaussian noise component
% lambda - Poisson parameter for the bursty signals (only relevant for
% simType = 3)
[signals, backgrounds] = generateSim(numTracks,numel(t),simType,expConst,noiseStd,2);

% Signal plotting controls
% To plot a few sample signals you can change these variables
% plotMin is the first signal plotted
% plotMax is the last signal plotted
% Example: to plot the first 5 signals change plotMax to 5
plotMin = 0;
plotMax = 5;


% Initialize cell arrays to store analysis results in case of comparison of
% stingency parameters
ampPosTotal = cell(length(stringencyParameter), 1);
ampNegTotal = cell(length(stringencyParameter), 1);
ampVarTotal = cell(length(stringencyParameter), 1);

eventsPerSignalTotal = cell(length(stringencyParameter), 1);       % Total events per signal

% Loop through each match probability
for SP = 1:length(stringencyParameter)
    % Perform analysis
    [timeBetweenPos, timeBetweenNeg, ampPos, ampNeg, ampVar, eventsPerSignal] = ...
        analysis(signals, backgrounds, timeVecTotal, ...
        movingAverageWindowSpan, stringencyParameter(SP), plotMin, plotMax);

    % Remove the -Inf placeholder used to initialize the variable
    ampPos = ampPos(2:end);
    ampNeg = ampNeg(2:end);
    ampVar = ampVar(2:end);
    
    timeBetweenPos = timeBetweenPos(2:end);
    timeBetweenNeg = timeBetweenNeg(2:end);
    
    % Store various analysis results if comparing between different
    % stringency parameters
    ampPosTotal{SP} = ampPos;
    ampVarTotal{SP} = ampVar;
    ampNegTotal{SP} = ampNeg;

    eventsPerSignalTotal{SP} = eventsPerSignal; 

    % Plot histograms and statistical analysis results
    figure;
    yyaxis left
    histogram(ampVar(abs(ampVar) < 250), 25);
    ylabel('Counts');
    yyaxis right
    var_gauss = pdf('Normal', -250:0.5:250, 0, std(ampVar(abs(ampVar) < 100)));
    plot(-250:0.5:250, var_gauss);
    ylabel('PDF');
    xlabel('Intensity [A.U]');
    title('Quiescent Amplitudes');


    % Plot histogram of all amplitudes with 3 Gaussians overlaid corresponding
    % to Negative, Zero, Positive

    figure; 
    histogram(ampPos, 'FaceColor', 'g', 'BinWidth', 5, 'Normalization', 'probability'); hold on;
    histogram(ampVar, 'FaceColor', 'b', 'BinWidth', 5, 'Normalization', 'probability');
    histogram(ampNeg, 'FaceColor', 'r', 'BinWidth', 5, 'Normalization', 'probability');
    xlim([-600, 600]); 
    ylabel('Frequency');
    title('3 Amplitude phases');
    xlabel('Amplitudes [A.U]');
    
    
    lambda = 1:10;
    K0Vector = 1:500;
    % Poisson analysis on positive amplitudes
    modifiedPoissonFit(ampPos,lambda,K0Vector,10,'Positive amplitudes');
    % Poisson analysis on negative amplitudes
    modifiedPoissonFit(abs(ampNeg),lambda,K0Vector,10,'Negative amplitudes');

end


