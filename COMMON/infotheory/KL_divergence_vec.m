function [Dkl]=KL_divergence_vec(distr1,distr2)
% changed name back from GEN_divergence_vec 29-11-2012
%[Dkl]=KL_divergence_vec(distr1,distr2)
%
%vector variant of KL_divergence
%calculates the Kullback-Leibler divergence of distr1 from distr 2.
%rows are considered as distributions
%gives a vector of one Dkl for each row
%Note: Dkl(A||B) is not Dkl(B||A), not a true distance, not symmetric
%extra preprocessing trick for skill scores:
%if distributions contain only one probability per row, it is assumed that that the
%distr is [1-p, p]
%
%note: distr1 and distr2 should be same size and distr2 should not contain
%zeros at places where distr1 is nonzero
%no further input-checking, use with caution
[nrows,ncols]=size(distr1);

%trick to convert prob of occurence of binary event into distribution
if ncols==1
    %assumes probabilities are between 0 and 1.
    p1=[1-distr1 distr1];
    p2=[1-distr2 distr2];
else
    %scale to probabilty distributions, in case histograms are used as inputs.
    p1 = distr1./repmat(sum(distr1,2),1,ncols);
    p2 = distr2./repmat(sum(distr2,2),1,ncols);
end




%calculate its entropy
%H=-sum((px.*log2(px)))
%note log2 is used to get info in unit "bits"


% %Dkl(size(p1))=0;
% %more robust version (including def 0*log2(0) = 0 )
% dDkl=0;
% for i=1:size(p1)
%     if p1(i)>0 
%         dDkl(i) = ( p1(i) * log2( p1(i) / p2(i) )  );
%     else
%         dDkl(i) = 0;
%     end
%     
% end
% 
% %Dkl is the Kullback-Leibler divergence
% Dkl = sum(dDkl);



dDkl= ( p1.* log2( p1 ./ p2 )  );
dDkl(p1==0)=0;

%Dkl is the Kullback-Leibler divergence
Dkl = sum(dDkl,2);
