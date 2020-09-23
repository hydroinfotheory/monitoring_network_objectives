function [H,Dkl]=histogram_entropy(histogram)
%calculates the shannon entropy of a histogram by scaling it to a discrete
%probability distribution. The kullback leibler distance from a uniform
%distribution is also calculated. The input is a vector with the number of
%elements per bin, that is the output of the "hist" function



%scale histogram to get discrete probability distribution
normalising_factor=sum(histogram);
%protection against empty histograms divide by zero
if normalising_factor==0
    H= 0;
    return
else
    px = histogram(:)/normalising_factor;
end

%calculate its entropy
%H=-sum((px.*log2(px)))
%note log2 is used to get info in unit "bits"


dH(size(px))=0;
%more robust version (including def 0*log2(0) = 0 )
H=0;
for i=1:size(px)
    if px(i)>0 
        dH(i) = -( px(i) * log2( px(i) )  );
    else
        dH(i) = 0;
    end
    
end

%H is Shannon Entropy (absolute, so independent of reference, but dependent
%on discretisation and does not converge if discr. approaches continuous
H = sum(dH); 

%Dkl is kullback-leibler divergence (independent on discretisation (lim to continuous converges,
%but is dependent on range of possible values (reference) because of log N  
%divergence is calculated relative to unform distribution over all discrete values
Dkl = log2(size(px,1))-H;

