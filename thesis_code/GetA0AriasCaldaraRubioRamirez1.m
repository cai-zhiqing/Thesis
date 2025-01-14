function [AA00,A0IN,irfs]=GetA0AriasCaldaraRubioRamirez1(bb,var,HorizonForSignRestrictions,LagOrder);
AA00=-9999;
A0IN=-9999;
irfs=-9999;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PositionRealInvestment=1;
PositionRealConsumption=2;
PositionRealImport=3;
PositionRealExport=4;
PositionLogM2=5;
PositionSpread=6;
PositionNEER=7;
PositionRealGDP=8;
PositionPrices=9;
PositionRepo7Day=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=size(var,1);
MaximumNumberOfRandomRotationMatrices=2000;
PositionMonetaryShock=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[VV,DD]=eig(var);
A0Start=VV*DD.^0.5;
InvA0StartPrime=(MyInverse(A0Start))';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumberOfZeroImpacts=7;
Z=zeros(NumberOfZeroImpacts,N);
Z(1:NumberOfZeroImpacts,1:NumberOfZeroImpacts)=eye(NumberOfZeroImpacts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IndicesMonetaryPolicyRule=[PositionRealGDP PositionPrices PositionRepo7Day]';
SignsMonetaryPolicyRule=[-1  -1  1]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
AA00=zeros(N,N,MaximumNumberOfRandomRotationMatrices);
A0IN=zeros(MaximumNumberOfRandomRotationMatrices,1);
%
zz=1;
while zz<=MaximumNumberOfRandomRotationMatrices
    Q=[];
    jj=1;
    while jj<=N
        if jj==PositionMonetaryShock
            RjofA0ApPlus=[Z*InvA0StartPrime; Q'];
        else
            RjofA0ApPlus=[Q'];
        end
        NN=null(RjofA0ApPlus);
        xj=randn(N,1);
        qj=NN*((NN'*xj)/norm(NN'*xj,2));
        Q=[Q qj];
        jj=jj+1;
    end
    %
    InvA0Prime=InvA0StartPrime*Q;
    %
    SignMonetaryPolicyRule=sign(InvA0Prime(IndicesMonetaryPolicyRule,PositionMonetaryShock));
    %
    if any(SignMonetaryPolicyRule-SignsMonetaryPolicyRule)==0
        IND=1;
    elseif any(-SignMonetaryPolicyRule-SignsMonetaryPolicyRule)==0
        InvA0Prime=-InvA0Prime;
        IND=1;
    else
        IND=0;
    end
    %
    if IND==1
        A0=MyInverse(InvA0Prime');
        %
        irfs=GetIRFs(bb,A0,N,LagOrder,HorizonForSignRestrictions,1);
        if mean(irfs(:,PositionRepo7Day)>0)==1 & mean(irfs(:,PositionSpread)<0)==1 & mean(irfs(:,PositionLogM2)<0)==1 & mean(irfs(:,PositionNEER)>0)==1
            AA00=A0;
            A0IN=1;
            break
        else
        end
    end
    zz=zz+1;
end



% IRFs(:,:,xx)=GetIRFs(BBBB(:,:,hh),AA00(:,:,hh),N,LagOrder,Horizon,1);
