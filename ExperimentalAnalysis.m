close all;   % Close all open figure windows
clear;       % Clear workspace

% Algorithm Parameters
movingAverageWindowSpan = 13;             % Moving average window span for signal denoising
stringencyParameter = 14;                 % Parameter for burst events detection
timeVecTotal = linspace(0, 60, 360);      % Time steps vector

% Load signal and background data
% To recreate the figures in the paper please uncomment any of the
% following signals,backgrounds pairs:

% RNA-protein granules formed with 1:1 protein to RNA ratio
signals = struct2cell(load('./in_vitro_experiments/protein_RNA_1_1_data.mat'));
backgrounds = struct2cell(load('./in_vitro_experiments/protein_RNA_1_1_background.mat'));

% RNA-protein granules formed with 10:1 protein to RNA ratio
% signals = struct2cell(load('./in_vitro_experiments/protein_RNA_10_1_data.mat'));
% backgrounds = struct2cell(load('./in_vitro_experiments/protein_RNA_10_1_background.mat'));

% RNA-protein granules formed with 100:1 protein to RNA ratio
% signals = struct2cell(load('./in_vitro_experiments/protein_RNA_100_1_data.mat'));
% backgrounds = struct2cell(load('./in_vitro_experiments/protein_RNA_100_1_background.mat'));

% In vivo granules formed with Qb-5x
% signals = struct2cell(load('./in_vivo_experiments/5Qb_data.mat'));
% backgrounds = struct2cell(load('./in_vivo_experiments/5Qb_background.mat'));

% In vivo granules formed with Qb-10x
% signals = struct2cell(load('./in_vivo_experiments/10Qb_data.mat'));
% backgrounds = struct2cell(load('./in_vivo_experiments/10Qb_background.mat'));


% Signal plotting controls
% To plot a few sample signals you can change these variables
% plotMin is the first signal plotted
% plotMax is the last signal plotted
% Example: to plot the first 5 signals change plotMax to 5
plotMin = 0;
plotMax = 0;


signals = signals{1};
backgrounds = backgrounds{1};
signals = signals(2:end);
backgrounds = backgrounds(2:end);


% Perform analysis
[timeBetweenPos, timeBetweenNeg, ampPos, ampNeg, ampVar, eventsPerSignal] = ...
        analysis(signals, backgrounds, timeVecTotal, ...
        movingAverageWindowSpan, stringencyParameter, plotMin, plotMax);

% Remove the -Inf placeholder used to initialize the variable
ampPos = ampPos(2:end);
ampNeg = ampNeg(2:end);
ampVar = ampVar(2:end);

timeBetweenPos = timeBetweenPos(2:end);
timeBetweenNeg = timeBetweenNeg(2:end);    

figure;
histogram(ampVar(abs(ampVar)<1500),50);
ylabel('Counts');
xlabel('Intesity [A.U]');
title('Quiescent Amplitudes');


% Plot histogram of all amplitudes 
figure; 
histogram(ampVar(abs(ampVar)<3000),'FaceColor','b','BinWidth',100);
hold on;
histogram(ampPos(ampPos<3000),'FaceColor','g','BinWidth',100); 
histogram(ampNeg(abs(ampNeg)<3000),'FaceColor','r','BinWidth',100);
ylabel('Counts'); xlim([-3000,3000])
title('Segment amplitudes','FontSize',12);
xlabel('Amplitudes [A.U]','FontSize',12);

% Plot boxplots of durations between burst events 
figure;
plot_var = [timeBetweenPos;timeBetweenNeg];
group_var = zeros(length(plot_var),1);
group_var(1:length(timeBetweenPos))=2;
group_var(length(timeBetweenPos)+1:end)=3;
boxplot(plot_var,group_var);
xticklabels({'Positive','Negative'});
ylabel('Time [min]');
title('Time between burst events');


lambda = 1:10;
K0Vector = 1:500;
% Poisson analysis on positive amplitudes
modifiedPoissonFit(rmoutliers(ampPos,"percentiles",[0,90]),lambda,K0Vector,10,'Positive amplitudes');
% Poisson analysis on negative amplitudes
modifiedPoissonFit(rmoutliers(abs(ampNeg),"percentiles",[0,90]),lambda,K0Vector,10,'Negative amplitudes');
