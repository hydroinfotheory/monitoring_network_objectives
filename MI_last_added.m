function [MIlast]=MI_last_added(data,col_order)
%  [MIlast]=MI_last_added(data,col_order)
%
%calculates the mutual information of the last station in col order with
%the combined information in the other stations in col order
% input
%   data :  quantized data, N rows (time steps), M columns (variables)
%   col_order: order in which the columns of the data are processed 
% output:
%   MIlast: mutual information of the time series is the last column with
%       the combined time series in the other columns.

k=0;
Hmarg=zeros(1,length(col_order));
JEk=Hmarg;
scum=data(:,col_order(1));  % re-organise data in column order of added stations 
for c=col_order             % loop over column is order of addition
    k=k+1;                  % counter
    s1=data(:,c);           % new added time series
    scum=merge_series(scum,s1);                     %add new series to combined multivariate merged series.
    Hmarg(k)=histogram_entropy(hist(s1,max(s1)));   %Calculate (marginal) entropy of added series 
    JEk(k)=histogram_entropy(hist(scum,max(scum))); %Calculate joint entropy of combined selected series
         
end
%mutual information is calculated by comparing the marginal entropy of the last station
%with the increase in joint entropy caused by adding it to the set
MIlast=Hmarg(end)-(JEk(end)-JEk(end-1));

