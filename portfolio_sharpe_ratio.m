function [results] = portfolio_sharpe_ratio(csv_dir, varargin)
% RESULTS = PORTFOLIO_SHARPE_RATIO(CSV_DIR, varargin)
% 
% *Info:
% This function optimizes portfolio weights based on Sharpe ratio,
% average total return, and volatility (standard deviation of returns). The
% basic idea is to provide a directory as input 'csv_dir'. This folder
% should contain .csv files of historical data for each ticker desired to 
% comprise part of a portfolio. The function will then return an optimized 
% weighting scheme based on high Sharpe ratio, high total return, and low 
% volatility (you can provide a custom weighting for each criterion). It 
% will also output historical performance data for the portfolio, an 
% alternative equally-weighted portfolio, and each individual ticker. 
% Furthermore, you can choose to plot a matrix showing correlations among 
% the individual assets, as well as a 3D scatter plot showing where the 
% optimized portfolio falls among other randomly generated portfolios on 
% the dimensions of Sharpe ratio, mean annual return, and volatility (SD of
% returns). The function was designed using data from Yahoo Finance
% (https://finance.yahoo.com/) but should work with other data sources
% provided the formatting is similar. Try optimizing the sample portfolio 
% included (with data from 2006-2018) to gain a better understanding of use
% cases.
% 
% *Cautionary note: your results may be skewed if your data does not go back
% sufficiently far and includes only one portion of a market cycle (e.g.,
% all bull market, no recessions). To gain insight into performance over
% the whole market cycle, try to include data going back to 2008 or 2000, 
% if not longer. When this is not possible, note that assets with high
% volatility and high annual returns during a bull market will often be the
% same assets that sustain the largest losses during economic downturns.
% In the case of major recessions these losses can sometimes exceed -40% 
% in a calendar year.
% 
% *Disclaimer: This open-source research tool is not intended to provide 
%   investment advice. It is intended only for informational purposes, and 
%   the user is not recommended to use the tool to make actual investment 
%   decisions. Seek a duly licensed professional for investment advice.
% 
% Inputs:
% 
%   'csv_dir',        a full path which contains .csv files of financial 
%                     data; files should be titled using the ETF/stock 
%                     ticker name
% 
% Optional Name-Value Pair Arguments:
% 
%   'nRandom',        "optimization" uses the simple method of iterating 
%                     over a large number of random portfolio weights and
%                     then selecting the best based on given criteria. The
%                     larger this number, the better the selection that
%                     can be made. (default:  10,000)
% 
%   'outcomeWeights', a 1x3 vector containing weights which
%                     correspond to optimization preferences for
%                     [large Sharpe ratio, large total return, small SD]
%                     (i.e., setting this parameter allows one to gear the
%                     portfolio optimization toward a more aggressive or
%                     conservative portfolio (default:  [1,1,1])
% 
%   'useGeoMeans',    boolean denoting whether to average across yearly
%                     returns and Sharpe ratios using geometric means 
%                     (true) or arithmetic means (false) (default: true)
% 
%   'minWeight',      numeric value (range: 0-1) specifying the minimum 
%                     weight that any ticker may receive (default: 0)
% 
%   'maxWeight',      numeric value (range: 0-1) specifying the maximum 
%                     weight that any ticker may receive (default: 1)
% 
%   'limitTickers',   numeric value (range: 0-nTickers) specifying the
%                     number of tickers to include in the optimized
%                     portfolio (which can be a subset of the number
%                     possible) (default:  include all)
% 
%   'plotScatter',    boolean denoting whether to plot a 3D scatter plot
%                     showing the optimized portfolio (red sphere) among
%                     randomly permuted portfolios on the dimensions of
%                     Sharpe, return, and SD (default: true)
% 
%   'plotCorrs',      boolean denoting whether to plot an imagesc plot
%                     showing the correlations between the returns of each
%                     ticker across all years (default: true)
% 
% Outputs:
% 
%   'results',        a structure containing:
% 
%       'results.weights'
%                     a table showing tickers with optimized weights
% 
%       'results.stats'
%                     a table of average Sharpe, total return, and SD 
%                     across years for the optimized portfolio, an equally 
%                     weighted portfolio, and ranges across permutations 
%                     for each measure
% 
%       'results.sharpe'
%                     a table displaying Sharpe ratio by year for both 
%                     the optimized portfolio, an equally weighted 
%                     portfolio, and for each individual ticker; tickers
%                     are sorted in descending order of average Sharpe
% 
%       'results.returns'
%                     the same as (.sharpe) but for total return; tickers
%                     are sorted in descending order of average returns
% 
%       'results.SD'
%                     the same as (.sharpe) but for SD; tickers are 
%                     sorted in ascending order of average SD
% 
% Example function call:
% 
%   results = portfolio_sharpe_ratio(folderPath,...
%       'nRandom',10000,'outcomeWeights',[.3,.3,.4],'useGeoMeans',1,...
%       'minWeight',.01,'maxWeight',.9,'limitTickers',10,'plotScatter',0,...
%       'plotCorrs',0);
% 
%   -this will use 10,000 randomly permuted portfolio weights, weight low
%   SD as a more important evaluation criterion than Sharpe or return, will
%   use geometric means, will not weight any portfolio asset less than 1%,
%   will not weight any portfolio asset more than 90%, will limit the
%   number of assets in the optimized portfolio to 10, and will not plot 
%   the 3D scatter plot or corrlation matrix.
% 
% Author:  Elliot Layden, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hWait1 = waitbar(.5,'Please wait...');

% Identify Function Path and Add Helper Scripts:
script_fullpath = mfilename('fullpath'); 
[script_path,~,~] = fileparts(script_fullpath); 
addpath(genpath(script_path))

% Get Inputs:
inputs = varargin;
parsed = struct('nRandom', 10000, 'outcomeWeights', ones(1,3), ...
    'useGeoMeans', true, 'minWeight', 0, 'maxWeight', 1, ...
    'limitTickers', 0, 'plotScatter', true, 'plotCorrs', true);
input_types = {{'numeric'},{'vector'},{'logical';'numeric'},...
    {'numeric'},{'numeric'},{'numeric'},{'logical';'numeric'},...
    {'logical';'numeric'}};  % column of internal cells == OR, row of internal cells == AND
parsed = getInputs(inputs, parsed, input_types);
nRandom = parsed.nRandom; outcomeWeights = parsed.outcomeWeights;
useGeoMeans = parsed.useGeoMeans; limitTickers = parsed.limitTickers;
minWeight = parsed.minWeight; maxWeight = parsed.maxWeight;
plotScatter = parsed.plotScatter; plotCorrs = parsed.plotCorrs;
if nRandom <=0; nRandom = 10000; end;
if length(outcomeWeights)~=3; outcomeWeights = ones(1,3); end

% Initialize Results:
results = struct('weights',[],'stats',[],'sharpe',[],'returns',[],'SD',[]);

% Extract portfolio data:
warning('off','all')
data = extract_portfolio(csv_dir);
warning('on','all')

% Get risk-free rates (10-year treasury rates):
warning('off','all')
riskFreeData = readtable(fullfile(script_path,'util','TNX.csv'));
warning('on','all')
riskFreeDates = year(riskFreeData.Date); yearsRiskFree = unique(riskFreeDates)';
nYear = length(yearsRiskFree); riskFree = zeros(1,nYear);
for i = 1:nYear
    vec = str2double(riskFreeData.Close(riskFreeDates==yearsRiskFree(i)));
    vec(isnan(vec)) = [];
    if useGeoMeans
        riskFree(i) = geomean(vec+1); riskFree(i) = riskFree(i) - 1;
    else
        riskFree(i) = mean(vec);
    end
end
riskFree(yearsRiskFree<year(data(1).date(1))) = []; riskFree = riskFree*.01;

% Extract all tickers:
nRows = size(data(1).return,1); nTickers = length(data);
dates = year(data(1).date); returns = zeros(nRows,nTickers); raw = returns;
tickers = cell(1, nTickers);
for i = 1:nTickers
   returns(:,i) = data(i).return; raw(:,i) = data(i).raw; tickers{i} = data(i).ticker;
end

% Input Checking:
if limitTickers >= nTickers; limitTickers = 0; end
if minWeight < 0; minWeight = 0; end
if maxWeight > 1; maxWeight = 1; end

% Permuted weights:
years = unique(dates)'; nYears = length(years); sharpe = zeros(1,nYears); 
yearlyReturns = sharpe; SD = sharpe; 
if limitTickers > 0
    weights = randfixedsum(limitTickers,nRandom,1,minWeight,maxWeight)';
    useTickers = zeros(limitTickers, nRandom);
else
    weights = randfixedsum(nTickers,nRandom,1,minWeight,maxWeight)';
    useTickers = repmat((1:nTickers)', 1, nRandom);
end
meanSharpe = zeros(1,nRandom); meanReturn = zeros(1,nRandom); meanSD = zeros(1,nRandom);
delete(hWait1); hWait = waitbar(0,'Processing...');
for perm = 1:nRandom
    waitbar(perm/nRandom, hWait);
    if limitTickers > 0
        useTickers(:,perm) = randperm(nTickers, limitTickers);
    end
    for i = 1:nYears
        [sharpe(i), yearlyReturns(i), SD(i)] = calculate_sharpe(raw(dates==years(i),useTickers(:,perm)), ...
            returns(dates==years(i),useTickers(:,perm)), weights(perm,:), riskFree(i));
    end
    if useGeoMeans
        meanSharpe(perm) = geomean(sharpe +10); 
        meanSharpe(perm) = meanSharpe(perm) -10;
        meanReturn(perm) = geomean(yearlyReturns + 10);
        meanReturn(perm) = meanReturn(perm) -10;
        meanSD(perm) = geomean(SD);
    else
        meanSharpe(perm) = mean(sharpe); 
        meanReturn(perm) = mean(yearlyReturns);
        meanSD(perm) = mean(SD);
    end
end; delete(hWait);

% Calculate minimum Euclidean distance from weighted best measures:
best = repmat([max(meanSharpe), max(meanReturn), min(meanSD)],nRandom,1);
check = [meanSharpe',meanReturn',meanSD'];
dist = sqrt(sum(repmat(outcomeWeights,nRandom,1).*(check - best).^2,2));
[~, ix] = min(dist); bestWeights = weights(ix,:); 
[bestWeights_sort, ixWeights] = sort(bestWeights,'descend');

% Add to results struct:
results.weights = cell2table(num2cell(round(bestWeights_sort*100,3)));
tickers1 = tickers(useTickers(:,ix)); tickers1 = tickers1(ixWeights);
results.weights.Properties.VariableNames = tickers1;
results.weights.Properties.RowNames = {'% weight:'};

% Calculate yearly for optimized weights:
sharpe_opt = zeros(1,nYears); yearlyReturns_opt = sharpe_opt; SD_opt = sharpe_opt;
for i = 1:nYears
    [sharpe_opt(i), yearlyReturns_opt(i), SD_opt(i)] = calculate_sharpe(raw(dates==years(i),useTickers(:,ix)),...
        returns(dates==years(i),useTickers(:,ix)), bestWeights, riskFree(i));
end    
% Calculate for equal weighting:
equalWeights = ones(1,nTickers)./nTickers; 
sharpe_eq = zeros(1,nYears); yearlyReturns_eq = sharpe_eq; SD_eq = sharpe_eq;
for i = 1:nYears
    [sharpe_eq(i), yearlyReturns_eq(i), SD_eq(i)] = calculate_sharpe(raw(dates==years(i),:), ...
        returns(dates==years(i),:), equalWeights, riskFree(i));
end
% Calculate stats for individual tickers:
sharpe_ind = zeros(nTickers, nYears); yearlyReturns_ind = sharpe_ind; SD_ind = sharpe_ind;
for i = 1:nTickers
    for j = 1:nYears
        [sharpe_ind(i,j), yearlyReturns_ind(i,j), SD_ind(i,j)] = calculate_sharpe(raw(dates==years(j),i), ...
            returns(dates==years(j),i), 1, riskFree(j));
    end
end

% Save table output 'stats':
if useGeoMeans
    results.stats = cell2table(num2cell(round([geomean(sharpe_opt +10), ...
        geomean(yearlyReturns_opt + 10), geomean(SD_opt);...
        geomean(sharpe_eq +10), geomean(yearlyReturns_eq +10), geomean(SD_eq);...
        max(meanSharpe),max(meanReturn),min(meanSD);...
        min(meanSharpe),min(meanReturn),max(meanSD)],3)));
    results.stats.Var1(1:2) = results.stats.Var1(1:2) - 10;
    results.stats.Var2(1:2) = results.stats.Var2(1:2) - 10;
else
    results.stats = cell2table(num2cell(round([mean(sharpe_opt), ...
        mean(yearlyReturns_opt), mean(SD_opt);...
        mean(sharpe_eq), mean(yearlyReturns_eq), mean(SD_eq);...
        max(meanSharpe),max(meanReturn),min(meanSD);...
        min(meanSharpe),min(meanReturn),max(meanSD)],3)));
end
results.stats.Properties.VariableNames = {'Sharpe','Return','SD'};
results.stats.Properties.RowNames = {'Optimized','EqualWeights','Best','Worst'};

yearLabels = sprintfc('%d',years);
for i = 1:length(yearLabels)
    yearLabels{i} = ['y',yearLabels{i}];
end
rowLabels = [{'Optimized','EqualWeights'},tickers];

% Save table output 'sharpe':
sharpeMat = [sharpe_opt; sharpe_eq; sharpe_ind];
if useGeoMeans
    [~, sharpeMat_ix] = sort(geomean(sharpeMat + 10,2),'descend');
else
    [~, sharpeMat_ix] = sort(mean(sharpeMat,2),'descend');
end
results.sharpe = cell2table(num2cell(round(sharpeMat(sharpeMat_ix,:),3)));
results.sharpe.Properties.VariableNames = yearLabels;
results.sharpe.Properties.RowNames = rowLabels(sharpeMat_ix);

% Save table output 'returns':
yearlyReturnsMat = [yearlyReturns_opt; yearlyReturns_eq; yearlyReturns_ind];
if useGeoMeans
    [~, yearlyReturns_ix] = sort(geomean(yearlyReturnsMat + 10,2),'descend');
else
    [~, yearlyReturns_ix] = sort(mean(yearlyReturnsMat,2),'descend');
end
results.returns = cell2table(num2cell(round(yearlyReturnsMat(yearlyReturns_ix,:),3)));
results.returns.Properties.VariableNames = yearLabels;
results.returns.Properties.RowNames = rowLabels(yearlyReturns_ix);

% Save table output 'SD':
SDMat = [SD_opt; SD_eq; SD_ind];
if useGeoMeans
    [~, SD_ix] = sort(geomean(SDMat,2),'ascend');
else
    [~, SD_ix] = sort(mean(SDMat,2),'ascend');
end
results.SD = cell2table(num2cell(round(SDMat(SD_ix,:),3)));
results.SD.Properties.VariableNames = yearLabels;
results.SD.Properties.RowNames = rowLabels(SD_ix);

% Plot Sharpe x Return x SD
if plotScatter
    if nRandom > 10000
        sparseIx = randperm(nRandom,10000);
    else sparseIx = 1:nRandom;
    end
    fig_title = 'Sharpe x Return x SD';
    figure('Name',fig_title,'units','norm','Color',[1,1,1],'Position',...
        [.061,.116,.861,.806],'NumberTitle','off','MenuBar','none'); 
    set(gca,'FontSize',24)
    scatter3(meanSharpe(sparseIx),meanReturn(sparseIx),meanSD(sparseIx),'.'); hold on;
    scatter3(meanSharpe(ix), meanReturn(ix), meanSD(ix),30,'o','fill','r')
    xlabel('Sharpe'); ylabel('Return'); zlabel('SD'); rotate3d
end

% Plot ticker correlation matrix:
if plotCorrs
    nTickers = size(returns,2); r = corr(returns); r(eye(nTickers)==1) = 0;
    fig_title = 'Ticker Correlation Matrix';
    figure('Name',fig_title,'units','norm','Color',[1,1,1],'Position',...
        [.236,.106,.554,.806],'NumberTitle','off','MenuBar','none'); 
    imagesc(r); 
    set(gca,'XTick',1:nTickers,'XTickLabels',tickers,'YTick',1:nTickers,...
        'YTickLabels',tickers,'XTickLabelRotation',90,'FontSize',12)
    colormap('jet'); colorbar; 
end

end