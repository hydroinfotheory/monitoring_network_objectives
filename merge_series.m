function s12=merge_series(s1,s2)
%merges 2 discretized series into one while preserving information
% Then relabels series to have discrete values from 1 to the number of
% unique values
%
% s12=merge_series(s1,s2)
%
% input : s1, s2 two discretized timeseries (vectors of positive integers)
% output: s12 a new series that has unique values for every encountered
%           combination of values for s1 and s2
%
%
% The entropy of s12 is equal to the joint entropy of s1, s2
ms1=max(s1);
raws12=ms1*(s2-1)+s1;
aa=unique(raws12);
s12=0*raws12;

%relabeling
k=0;
for val=aa'
    k=k+1;
    s12(raws12==val)=k;
end

