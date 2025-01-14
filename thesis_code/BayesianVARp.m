function [B,V,U]=BayesianVARp(X,LagOrder,Intercept,NumberOfDraws,Stationary);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luca Benati
% December 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [B,V,U]=BayesianVARp(X,LagOrder,Intercept,NumberOfDraws,Stationary);
% This program estimates a Bayesian VAR(p) for K series, stored in the
% matrix X. It is based on the formulas found in Uhlig (Carnegie-Rochester,
% 1998), and Uhlig (JME, 2005).
%                                             Input of the program is:
%
% X             = a matrix containing the K series stacked in columns: first column,
%                 the first series, second column, the second series, etc.
% LagOrder      = the order of the VAR
% Intercept     = 'Y' if you want an intercept vector, 'N' if you don't want it
% NumberOfDraws = the number of draws from the posterior distribution
% Stationary    = 'Y' if you want to impose a stationarity constraint on the VAR, 'N' otherwise
%
%                                             Output of the program is:
%
% B             = the posterior distribution of the estimated VAR,
%                 organized as a 3-dimensional object with each 'slice' of the object
%                 being structured as [MU B(1) B(2) ... B(p)], or as [MU
%                 DeterministicTimeTrend B(1) B(2) ... B(p)] in case in
%                 which there is also a deterministic time trend
% V             = the posterior distribution of the VAR's variance-covariance
%                 matrix, organized as a 3-dimensional object
% U             = the posterior distribution of the VAR's residuals,
%                 organized as a 3-dimensional object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargin<6
    Stationary='N';
end
%
[T,K]=size(X);
%
Y=X(LagOrder+1:T,:)';
%
if Intercept=='Y'
    Z=ones(1,T-LagOrder);
else
    Z=[];
end
for kk=1:LagOrder
    Z=[Z; X(LagOrder+1-kk:T-kk,:)'];
end
% OLS estimator:
B=(Y*Z')*MyInverse(Z*Z'); % Expression (3.2.10) in Lutkepohl
% Residuals:
U=Y-B*Z;
% Variance-covariance matrix:
V=(U*U')/(T-K*LagOrder-1);  % Expression (3.2.19) in Lutkepohl
%
if NumberOfDraws==0
    Sigma=V;
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The posterior (see Uhlig, JME, 2005, p. 410):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
InverseV=MyInverse(V);
%
Sigma=zeros(K,K,NumberOfDraws);
xx=1;
while xx<=NumberOfDraws
    Sigma(:,:,xx)=MyInverse(WishartRandomDraw(InverseV/T,T));
    xx=xx+1;
end
%
vecBT=vec(B');
InverseNT=MyInverse(Z*Z');
%
B=zeros(size(B,1),size(B,2),NumberOfDraws);
U=zeros(size(U,1),size(U,2),NumberOfDraws);
xx=1;
while xx<=NumberOfDraws
    CovarianceMatrix=kron(Sigma(:,:,xx),InverseNT);
    vecB=MultivariateNormalRandomDraw(vecBT,CovarianceMatrix,1)';
    bb=reshape(vecB,size(B,2),size(B,1))';
    if Stationary=='N'
        B(:,:,xx)=bb;
        U(:,:,xx)=Y-bb*Z;
        xx=xx+1;
    else
        if Intercept=='Y'
            MaxAbsRoot=max(abs(varroots(LagOrder,K,bb(:,1:1+LagOrder*K))));
        else
            MaxAbsRoot=max(abs(varroots(LagOrder,K,[zeros(K,1) bb(:,1:LagOrder*K)])));
        end
        if abs(MaxAbsRoot-1)<(0.5/T)
            B(:,:,xx)=bb;
            U(:,:,xx)=Y-bb*Z;
            xx=xx+1;
        end
    end
end
%
V=Sigma;















