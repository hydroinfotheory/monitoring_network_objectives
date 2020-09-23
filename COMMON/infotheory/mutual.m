function [MI,JE,Hcond_x_y,Hcond_y_x] = mutual(joint_distribution)
% [MI,JE,Hcond_x_y,Hcond_y_x] = mutual(joint_distribution)
% calculates various information-theoretical properties from a joint pdf
% INPUTS: 
%       joint_distribution : a matrix of probabilities (or densities scaled
%       to sum to one)
%       the different rows represent x value intervals, the columns y value intervals 
% OUTPUTS:
%       MI          : Mutual information between x and y
%       JE          : Joint Entropy of x and y
%       Hcond_x_y   : Conditional entropy of x , given y
%       Hcond_y_x   : Conditional entropy of y , given x

jd=joint_distribution;
nr_x= size(jd,1); %rows
nr_y= size(jd,2); %cols

%rescale jd to make sum 1
total_sum=sum(sum(jd));
jd=jd./total_sum;

for i=1:nr_x
    marg_px(i)= sum(jd(i,:));
    H_y_given_x(i)=histogram_entropy(jd(i,:));
end
Hcond_y_x = (marg_px*H_y_given_x');
Hx =  histogram_entropy(marg_px);

for j=1:nr_y
    marg_py(j)=sum(jd(:,j));
    H_x_given_y(j)=histogram_entropy(jd(:,j));
end
Hcond_x_y = (marg_py*H_x_given_y');
Hy =  histogram_entropy(marg_py);

MI1=Hy-Hcond_y_x;
MI2=Hx-Hcond_x_y;

for i=1:nr_x
   total_vector( 1+((i-1)*nr_y):i*nr_y)=jd(i,:);
end


joint_entropy=histogram_entropy(total_vector);

%corrected MI
%MI_corr = MI + (( (nr_x*nr_y)-nr_x-nr_y+1 ) / (2*samples  ) )

%outputs
MI=MI1;
JE=joint_entropy;
Hcond_x_y;
Hcond_y_x;
