function [sharpe, totalReturn, portfolioSD] = calculate_sharpe(raw, returns, weights, riskFreeReturn)
    SD = std(returns) .* sqrt(size(returns,1)); % annualized SD's
    wSD = weights .* SD;
    portfolioSD = sqrt(wSD * corr(returns) * wSD');
    totalReturn = (raw(end,:) - raw(1,:)) ./ raw(1,:);
    totalReturn = sum(totalReturn.*weights);
    sharpe = (totalReturn - riskFreeReturn) / portfolioSD;
end