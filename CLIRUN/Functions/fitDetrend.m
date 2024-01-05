function [detrendedData, p] = fitDetrend(data, daterange)
% function fitDetrend: takes in GCMs Precip or Temp data and a daterange
% fits a fourth order polynomial and then removes the trend from the data
% input:
%   data: GCMs Precipitation or Temperature data
%   daterange: daterange for GCMs data
% output
%   detrendedData: data without the trend (not including y-intercept)

% get lengths of timeseries and number of samples
numsamp = size(data, 2);
T = length(daterange);
t = (1:T)';
X = [t.^4 t.^3 t.^2 t ones(T,1)];
p = zeros(numsamp, size(X, 2));

for i=1:numsamp
    % get time series
    timeSeries = data(daterange,i);
    
    % use polyfit to fit 4th order polynomial (Hawkins and Sutton 2009)
    p(i,:) = polyfit(t, timeSeries, 4);
    
    % detrend data
    b = p(i,:)';
    tH = X*b;
    y_int = repmat(p(i,end),T, 1);
    detrendedData(:,i) = timeSeries - tH + y_int;
    
    
end
end