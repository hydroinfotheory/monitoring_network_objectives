function [Dkl]=KL_divergence(distr1,distr2)
%calculates the Kullback-Leibler divergence of distr1 from distr 2. 
%Note: Dkl(A||B) is not Dkl(B||A), not a true distance, not symmetric

%scale to probabilty distributions, in case historgams are used as inputs.
p1 = distr1(:)/sum(distr1);
p2 = distr2(:)/sum(distr2);


%calculate its entropy
%H=-sum((px.*log2(px)))
%note log2 is used to get info in unit "bits"


%Dkl(size(p1))=0;
%more robust version (including def 0*log2(0) = 0 )
dDkl=0;
for i=1:size(p1)
    if p1(i)>0 
        dDkl(i) = ( p1(i) * log2( p1(i) / p2(i) )  );
    else
        dDkl(i) = 0;
    end
    
end

%Dkl is the Kullback-Leibler divergence
Dkl = sum(dDkl);