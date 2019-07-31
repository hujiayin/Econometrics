% Pre-Estimation Analysis
% 1. Load the raw data: daily exchange rate
Close
T = length(Close);
plot([0:T-1], Close)
set(gca, 'XTick', [1 252 503 756 1008 T])
set(gca, 'XTickLabel', {'Feb 2014' 'Feb 2015' 'Feb 2016' 'Feb 2017' ... 
    'Feb 2018' 'Feb 2019'})
ylabel('S&P 500 Index Prices')
title('S&P 500 Index Prices 2014/02/03 - 2019/02/01')

% 2. Convert the prices to a return series
SPdata_r = price2ret(Close);

plot(SPdata_r)
set(gca,'XTick',[1 251 502 755 1007 T-1]) 
set(gca,'XTickLabel',{'Feb 2014' 'Feb 2015' 'Feb 2016' 'Feb 2017' ... 
    'Feb 2018' 'Feb 2019'})
ylabel('Return')
title('S&P 500 Index Prices Daily Returns 2014/02/03 - 2019/02/01')

% 3. Check for correlation in the return series
autocorr(SPdata_r)
title('ACF with Bounds for Raw Return Series')

parcorr(SPdata_r)
title('PACF with Bounds for Raw Return Series')

% 4. Check for correlation in the squared returns
autocorr(SPdata_r.^2)
title('ACF of the Squared Returns')

% 5. Quantify the correlation.
% Ljung-Box-Pierce Q-test
[H,pValue,Stat,CriticalValue] = lbqtest(SPdata_r - mean(SPdata_r),'lag', [10 15 20]);
[H pValue Stat CriticalValue]

[H,pValue,Stat,CriticalValue] = lbqtest((SPdata_r-mean(SPdata_r)).^2,'lag', [10 15 20]);
[H pValue Stat CriticalValue]

% Engle¡¯s ARCH Test
[H,pValue,Stat,CriticalValue] = archtest(SPdata_r - mean(SPdata_r),'lag', [10 15 20]);
[H pValue Stat CriticalValue]

% Parameter Estimation
% 1. Estimate the Model Parameters
ToEstMdl = garch('Offset',NaN,'GARCHLags',1,'ARCHLags',1);
[EstMdl,EstParamCov,logL,info] = estimate(ToEstMdl, SPdata_r)

% Post-estimation Analysis
% 1. Compare the Residuals, Conditional Standard Deviations
% 2. Plot and Compare the Correlation of the Standardized Innovations
v = infer(EstMdl, SPdata_r);

figure
plot(v)
xlim([0,T])
title('Inferred Conditional Variances')

res = (SPdata_r-EstMdl.Offset)./sqrt(v);

figure
subplot(2,2,1)
plot(res)
xlim([0,T])
title('Standardized Residuals')

subplot(2,2,2)
histogram(res,10)

subplot(2,2,3)
autocorr(res)

subplot(2,2,4)
parcorr(res)

% 3. Quantify and Compare Correlation of the Standardized Innovations
[H, pValue,Stat,CriticalValue] = lbqtest(res,'lag', [10 15 20]);


% Garch(2,1)
% Parameter Estimation
% 1. Estimate the Model Parameters
ToEstMdl = garch('Offset',NaN,'GARCHLags',[1 2],'ARCHLags',1);
[EstMdl2,EstParamCov2,log2L,info2] = estimate(ToEstMdl, SPdata_r);


% Garch(1,2)
% Parameter Estimation
% 1. Estimate the Model Parameters
ToEstMdl = garch('Offset',NaN,'GARCHLags',1,'ARCHLags',[1 2]);
[EstMdl2,EstParamCov2,logL2,info2] = estimate(ToEstMdl, SPdata_r);



% Forecast
vF1 = forecast(EstMdl,1,'Y0',SPdata_r)

% Robustness
SPdata_r1 = SPdata_r(1:(T-1)/2);
SPdata_r2 = SPdata_r((T-1)/2 + 1: T-1);

ToEstMdl = garch('Offset',NaN,'GARCHLags',1,'ARCHLags', 1);
[EstMdl_1,EstParamCov_1,logL_1,info_1] = estimate(ToEstMdl, SPdata_r1);

ToEstMdl = garch('Offset',NaN,'GARCHLags',1,'ARCHLags', 1);
[EstMdl_2,EstParamCov_2,logL_2,info_2] = estimate(ToEstMdl, SPdata_r2);

