function FractionsOfVariance=GetFractionsOfVariance(B,A0,N,LagOrder,Horizon,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I: The variances conditional on all shocks:
[a,A,VAR]=Companion(B,A0*A0',N,LagOrder);
VariancesAllShocks=[];
VariancesOld=zeros(size(A));
kk=0;
while kk<=Horizon
    % If kk=0 it is on impact; if kk=1 it is one period ahead, etc.
    Akk=A^kk;
    VariancesNew=VariancesOld+Akk*VAR*Akk';
    VariancesAllShocks=[VariancesAllShocks; diag(VariancesNew(1:N,1:N))'];
    VariancesOld=VariancesNew;
    kk=kk+1;
end
%
FractionsOfVariance=zeros(Horizon+1,N,K);
for xx=1:K % Shock
    % II: The variances conditional on individual shocks:
    CovarianceMatrixStructuralShocks=zeros(N,N);
    CovarianceMatrixStructuralShocks(xx,xx)=1;
    [a,A,VAR]=Companion(B,A0*CovarianceMatrixStructuralShocks*A0',N,LagOrder);
    Variances=[];
    VariancesOld=zeros(size(A));
    kk=0;
    while kk<=Horizon
        Akk=A^kk;
        VariancesNew=VariancesOld+Akk*VAR*Akk';
        Variances=[Variances; diag(VariancesNew(1:N,1:N))'];
        VariancesOld=VariancesNew;
        kk=kk+1;
    end
    FractionsOfVariance(:,:,xx)=Variances./VariancesAllShocks;
end






