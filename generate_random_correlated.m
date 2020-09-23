function series=generate_random_correlated(corrmat,N,make_uniform)
%   series=generate_random_correlated(corrmat,N)
%generate correlated normally distributed series with specified correlation matrix as given by corrcoef
% series is N by length(corrmat) 
if nargin<3
    make_uniform=0;
end


%generate correlated random normal vars
U = chol(corrmat); %cholesky decomposition
series=randn(N,length(corrmat))*U; %generate vars
if make_uniform
    series=normcdf(series);
end

