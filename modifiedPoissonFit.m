% Function modifiedPoissonFit fits the amplitude data to modified Poisson
% functions, calcualtes the mean squared error between the observed data
% and the theoretical functions, and plots the results
function [pvec, MSE, QQ_error, KS_h, KS_dist, w_score] = modifiedPoissonFit(data,lambda,K0Vector,numBins,title_str)

% Present the data as 10 bins
[binCounts,binEdges] = histcounts(data,numBins,'Normalization','pdf');

if binEdges(1) < 0
    binEdges = fliplr(binEdges);
    binCounts = fliplr(binCounts);
end
nn_full = binCounts;
hist_x_axis = abs(mean([binEdges(1:end-1);binEdges(2:end)]));
percent = 15;
GOF_threshold = prctile(hist_x_axis, percent); % The pecrentile as a threshold value for goodness-of-fit analysis (we are interested more about the right side of distribution)
% If we want no data to be truncated, GOF_threshold & percent should be set to 0
% GOF_threshold = 0;percent = 0;

% Initialize MSE matrix
MSE = zeros(length(lambda),length(K0Vector));
QQ_error = zeros(length(lambda),length(K0Vector)); % QQ plot RMSE
KS_h = zeros(length(lambda),length(K0Vector)); % Decision of Kolmogorov-Smirnov test
KS_dist = zeros(length(lambda),length(K0Vector)); % Kolmogorov-Smirnov statistic

w_MSE = 0.5; % Weight for MSE
w_QQ = 0.25;  % Weight for QQ RMSE
w_KS = 0.25;  % Weight for KS distance

x_axis = zeros(1,length(K0Vector));

for j = 1:length(lambda)
    for i=1:length(K0Vector)
        normalizer = K0Vector(i); % Current normalization factor
        x_axis(1,i) = normalizer;
        xx = 0:max(hist_x_axis)/normalizer+1;
        
        % Generate a theoretical Poisson distribution using the normalized
        % x axis and the current lambda from the lambda vector
        yy = poisspdf(xx,lambda(j));
        yq = interp1(xx,yy./max(yy),linspace(0,max(hist_x_axis)/normalizer,length(binCounts)));
        
        % Calcualte goodness of fit between the truncated data and the theoretical function
        yq_thres = yq(hist_x_axis > GOF_threshold);
        binCounts_thres = binCounts((hist_x_axis > GOF_threshold));
        nn_full_thres = nn_full((hist_x_axis > GOF_threshold));

        if all(~yq_thres) || all(~nn_full_thres) || all(~binCounts_thres)
            warning(['Truncation unsuccesful for lambda=',num2str(lambda(j)),' and k0=',num2str(normalizer),' - working on full data'])
            yq_thres = yq;
            binCounts_thres = binCounts;
            nn_full_thres = nn_full;
        end
        
        MSE(j,i) = goodnessOfFit(yq_thres(~isnan(yq_thres)),binCounts_thres(~isnan(yq_thres))./max(binCounts_thres),'MSE');
        n = numel(nn_full_thres);
        p = ((1:n) - 0.5) / n;
        q1 = quantile(nn_full_thres, p);
        q2 = quantile(yq_thres*max(binCounts_thres)./max(yq_thres), p);
        QQ_error(j,i) = sqrt(mean((q2 - q1).^2));
        [~,~,KS_dist(j,i)] = kstest2(yq_thres*max(binCounts_thres)./max(yq_thres), nn_full_thres, 0.01); % Two-sample Kolmogorov-Smirnov test
    end
end

% Ranking Normalization
M = numel(lambda)*numel(K0Vector);
rankKS = tiedrank(KS_dist(:)); rankKS = (rankKS-1)/(M-1); rankKS = reshape(rankKS, size(KS_dist));
rankMSE = tiedrank(MSE(:)); rankMSE = (rankMSE-1)/(M-1); rankMSE = reshape(rankMSE, size(MSE));
rankQQ = tiedrank(QQ_error(:)); rankQQ = (rankQQ-1)/(M-1); rankQQ = reshape(rankQQ, size(QQ_error));
w_score = (w_QQ * rankQQ) + (w_KS * rankKS) + (w_MSE * rankMSE); % Weighted scores

fprintf('\n\n%s weights:\n', title_str)
fprintf('  %-25s : %g\n', "Weight MSE", w_MSE);
fprintf('  %-25s : %g\n', "Weight QQ", w_QQ);
fprintf('  %-25s : %g\n', "Weight KS", w_KS);
fprintf('  %-25s : %g\n\n', "Truncation Percentile", percent);


% Plot a subset of the results (since the K0Vector can be hundreds of values).
norm_short = K0Vector(1:20:end);
lambda_short = lambda;
w_short = w_score(:,1:20:end);

% Plot WS heatmap
figure;
heatmap(norm_short,lambda_short,w_short,'CellLabelColor','none'); colorbar();
colormap(parula); caxis([0,1]);
xlabel('Model single molecule intensity');
ylabel('Poisson rate');
title(['Poisson fitting parameters', newline, title_str]);

% plot (lambda,k0) fits with relevant QQ plots and CDF plots
[minScore,K_idx] = min(w_score,[],2);
for i = 1:length(lambda)
    figure;
    s = subplot(1,5,1:3);
    pvec(i) = s;
    scatter(hist_x_axis,nn_full,'DisplayName',[title_str,' experimental data']);  
    box on;
    hold on;
    lambda_single = lambda(1,i);
    normalizer = x_axis(1,K_idx(i));
    xx = 0:max(hist_x_axis)/normalizer+1;
    yy = poisspdf(xx,lambda_single);
    yq = interp1(xx,yy./max(yy),linspace(0,max(hist_x_axis)/normalizer,length(binCounts)));
    yq_thres = yq(hist_x_axis > GOF_threshold);
    binCounts_thres = binCounts((hist_x_axis > GOF_threshold));
    nn_full_thres = nn_full((hist_x_axis > GOF_threshold));

    if all(~yq_thres) || all(~nn_full_thres) || all(~binCounts_thres)
        warning(['Truncation unsuccesful for lambda=',num2str(lambda(j)),' and k0=',num2str(normalizer),' - working on full data'])
        yq_thres = yq;
        binCounts_thres = binCounts;
        nn_full_thres = nn_full;
    end

    legend_txt = [' K0=',num2str(normalizer)];
    plot(hist_x_axis,yq*max(binCounts)./max(yq),'--r','LineWidth',1.5,'DisplayName',legend_txt);
    xlabel('Intensity [A.U]');
    ylabel('Probability');
    yl = get(s, "YLim");
    xl = get(s, "XLim");
    text(xl(2)*0.75, yl(2)*0.8, ['\lambda_{est}=', num2str(lambda_single),' K_0=', num2str(normalizer),newline,' WS=',num2str(round(minScore(i),2,'significant'))], 'FontSize', 14, 'FontWeight', 'bold');
    legend('show');

    subplot(1,5,4);
    qqplot(yq_thres*max(binCounts_thres)./max(yq_thres),nn_full_thres); box on;
    title('Truncated QQ plot');
    ylabel('Quantiles of sample data');
    xlabel('Quantiles of theoretical distribution');

    subplot(1,5,5);
    cdfplot(nn_full_thres); hold on;
    cdfplot(yq_thres*max(binCounts_thres)./max(yq_thres));
    title('Truncated CDF plot');
    legend('Experiment','Estimate');
end



end
