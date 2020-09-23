function [JE,TC,Hmarg]=merging_stats_fast(data,col_order)
%
% [JE,TC,Hmarg]=merging_stats_fast(data,col_order)
%
%  Calculates information theoretical statistics for a group of timeseries
%  inputs: data - N by M matrix of discretized data (positive integers) 
%                 the N rows are timesteps, the M columns are the
%                 different time series
%          col_order - a vector on integers between 1 and M, indicating
%          which columns are used and in which order. 
%  
%  
%  outputs: JE -    The joint entropy calculated for the entire set of
%                   selected timeseries  data(:,col_order)
%           TC -    The total correction for that same set
%           Hmarg - A vector of marginal entropies for each selected column
%

k=0;
Hmarg=zeros(1,length(col_order));
% scum= is the cumulative set of selected time series, which is grown one
% by one in the loop below by merging the new gauge with the existing set
scum=data(:,col_order(1));
%loop 
for c=col_order
    k=k+1;
    s1=data(:,c);  % s1 is newest selected time series
    scum=merge_series(scum,s1); %merge existing selected series with new one
    Hmarg(k)=histogram_entropy(hist(s1,max(s1))); %calculate marginal entropy for each selected column individually
end
JE=histogram_entropy(hist(scum,max(scum)));   %entropy of merged time series is joint entropy of all contained therein
sumH=sum(Hmarg);                              %sum of marginal entropies would be joint entropy if all are indepedent
TC=sumH-JE;                                   %total correlation (multivariate extension of mutual information) is the difference
