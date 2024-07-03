% Function modifiedPoissonFit fits the amplitude data to modified Poisson
% functions, calcualtes the mean squared error between the observed data
% and the theoretical functions, and plots the results
function modifiedPoissonFit(data,lambda,K0Vector,numBins,title_str)

% Present the data as 10 bins
[binCounts,binEdges] = histcounts(data,numBins,'Normalization','pdf');

if binEdges(1) < 0
    binEdges = fliplr(binEdges);
    binCounts = fliplr(binCounts);
end
nn_full = binCounts;
hist_x_axis = abs(mean([binEdges(1:end-1);binEdges(2:end)]));

% Initialize MSE matrix
MSE = zeros(length(lambda),length(K0Vector));
mse_x_axis = zeros(1,length(K0Vector));
for j = 1:length(lambda)
    for i=1:length(K0Vector)
        normalizer = K0Vector(i); % Current normalization factor
        mse_x_axis(1,i) = normalizer;
        xx = 0:max(hist_x_axis)/normalizer+1;
        % Generate a theoretical Poisson distribution using the normalized
        % x axis and the current lambda from the lambda vector
        yy = poisspdf(xx,lambda(j));
        yq = interp1(xx,yy./max(yy),linspace(0,max(hist_x_axis)/normalizer,length(binCounts)));
        % Calcualte goodness of fit between the data and the theoretical
        % function
        MSE(j,i) = goodnessOfFit(yq(~isnan(yq)),binCounts(~isnan(yq))./max(binCounts),'MSE');
    end
end

% Plot a subset of the results (since the K0Vector can be hundreds of
% values).
norm_short = K0Vector(1:20:end);
lambda_short = lambda;
MSE_short = MSE(:,1:20:end);
% Plot MSE heatmap
figure;
heatmap(norm_short,lambda_short,MSE_short,'CellLabelColor','none'); colorbar();
colormap(parula); caxis([0,5]);
xlabel('Model single molecule intensity');
ylabel('Poisson rate');
title('Poisson fitting parameters');
minMatrix = min(MSE(:));
[row,col] = find(MSE==minMatrix,1);
lambda_single = lambda(1,row);
normalizer = mse_x_axis(1,col);
xx = 0:max(hist_x_axis)/normalizer+1;
yy = poisspdf(xx,lambda_single);
yq = interp1(xx,yy./max(yy),linspace(0,max(hist_x_axis)/normalizer,length(binCounts)));


% plot (lambda,k0) fits with relevant QQ plots and CDF plots
[minMatrix,minInd] = min(MSE,[],2);

for i = 1:5
    figure;
    subplot(1,5,1:3);
    scatter(hist_x_axis,nn_full,'DisplayName',[title_str,' experimental data']);  
    box on;
    hold on;
    lambda_single = lambda(1,i);
    normalizer = mse_x_axis(1,minInd(i));
    xx = 0:max(hist_x_axis)/normalizer+1;
    yy = poisspdf(xx,lambda_single);
    yq = interp1(xx,yy./max(yy),linspace(0,max(hist_x_axis)/normalizer,length(binCounts)));
    legend_txt = [' K0=',num2str(normalizer),' MSE=',num2str(minMatrix(i))];
    plot(hist_x_axis,yq*max(binCounts)./max(yq),'--r','LineWidth',1.5,'DisplayName',legend_txt);
    xlabel('Intensity [A.U]');
    ylabel('Probability');
    legend('show');
    title('lambda=',num2str(lambda_single));
    
    subplot(1,5,4);
    qqplot(yq*max(binCounts)./max(yq),nn_full); box on;
    title('QQ plot');
    ylabel('Quantiles of sample data');
    xlabel('Quantiles of theoretical distribution');
    
    subplot(1,5,5);
    cdfplot(nn_full); hold on;
    cdfplot(yq*max(binCounts)./max(yq));
    title('CDF plot');
    legend('Experiment','Estimate');
end


end
