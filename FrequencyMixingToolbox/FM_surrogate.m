function [X_s, goodRand] = FM_surrogate(X,k)

covX = cov(X);
%covX = X'*X;%/size(X,1);
chol_covX = chol(nearestSPD(covX));
temp = mvnrnd([0 0],[1 1],k*size(X,2));
rand1_ = temp(:,1); rand2_ = temp(:,2);
rand1 = reshape(rand1_,[k size(X,2)]);
rand2 = reshape(rand2_,[k size(X,2)]);
goodRand = complex(rand1,rand2)/sqrt(2); %want mu = 0, sigma = 1
%goodRand = complex(rand1,rand2)./abs(complex(rand1,rand2)); %want mu = 0, sigma = 1
%temp = complex(randn(k,size(X,2)), randn(k,size(X,2)));

X_s = [goodRand*chol_covX];



end

