%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
path(path,'U:\EmpiricalMacro')
options = optimset;
options = optimset(options,'Display','off');
warning('off','all')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumberOfDraws=5000;
LagOrder=6;
HorizonForSignRestrictions=24;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=xlsread('ChinaMonthlyData.xlsx','MonthlyData','A2:M325');
[T,N]=size(X);
Time=X(:,1);
LogRealGDP=log(X(:,2));
LogNomInvestment=log(X(:,3));
LogNomConsumption=log(X(:,4));
LogM2=log(X(:,5));
LogNomImport=log(X(:,6));
LogNomExport=log(X(:,7));
Repo7Day=X(:,8);
DepositRate1Year=X(:,9);
LogGDPDeflator=log(X(:,11));
LogCPI=log(X(:,12));
LogNEER=log(X(:,13));
LogRealInvestment=LogNomInvestment-LogGDPDeflator;
LogRealConsumption=LogNomConsumption-LogGDPDeflator;
LogRealImport=LogNomImport-LogGDPDeflator;
LogRealExport=LogNomExport-LogGDPDeflator;
Spread=DepositRate1Year-Repo7Day;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=[LogRealInvestment LogRealConsumption LogRealImport LogRealExport LogM2 Spread LogNEER LogRealGDP LogGDPDeflator Repo7Day];
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
%
[T,N]=size(X);
Stationary='N';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Trend=0;
Intercept='Y';
Horizon=10*12;
T=length(Time);
NN=NumberOfDraws;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[BBBB,SSSS,UUUU]=BayesianVARp(X,LagOrder,Intercept,NumberOfDraws,Stationary);
if Intercept=='Y'
    MUMU=BBBB(:,1,:);
    BBBB=BBBB(:,2:size(BBBB,2),:);
else
    MUMU=zeros(N,1,NN);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DFILE=['ThesisImposingRestriction'];
varname(1,:)=['BBBB'];
varname(2,:)=['MUMU'];
varname(3,:)=['SSSS'];
varname(4,:)=['UUUU'];
varname(5,:)=['AA00'];
varname(6,:)=['A0IN'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AA00=zeros(N,N,NumberOfDraws);
A0IN=zeros(NumberOfDraws,1);
hh=1;
while hh<=NumberOfDraws
    [A0,IND,irfs]=GetA0AriasCaldaraRubioRamirez1(BBBB(:,:,hh),SSSS(:,:,hh),HorizonForSignRestrictions,LagOrder);
    if IND==1
        AA00(:,:,hh)=A0;
        A0IN(hh)=IND;
        sum(A0IN)
    end
    hh=hh+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NNN=sum(sum(A0IN)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                Impulse-response functions
IRFs=zeros(Horizon+1,N,NNN);
%
xx=0;
hh=1;
while hh<=NumberOfDraws
    if A0IN(hh)==1
        xx=xx+1;
        IRFs(:,:,xx)=GetIRFs(BBBB(:,:,hh),AA00(:,:,hh),N,LagOrder,Horizon,1);
        IRFs(:,:,xx)=(IRFs(:,:,xx)/median(IRFs(1,PositionRepo7Day,xx),3))*0.25;
    end
    hh=hh+1;
end
SortedIRFs=sort(IRFs,3);
HOR=(0:1:Horizon)';
%
Percentiles=fix(NNN*[0.5 0.16 0.84 0.05 0.95]');
%
PercentilesIRFs=SortedIRFs(:,:,Percentiles);
%
for xx=1:N
    Perc=squeeze(PercentilesIRFs(:,xx,:));
    figure(1)
    subplot(2,N,xx)
    plot(HOR,zeros(size(HOR)),'b:',HOR,Perc(:,1),'k',HOR,Perc(:,2:3),'r','LineWidth',2)
    xlim([0 Horizon])
    if xx==PositionRealInvestment
        title('Investment')
        ylabel('IRFs')
    elseif xx==PositionRealConsumption
        title('Consumption')
    elseif xx==PositionRealImport
        title('Import')
    elseif xx==PositionRealExport
        title('Export')
    elseif xx==PositionLogM2
        title('M2')
    elseif xx==PositionSpread
        title('Spread')
    elseif xx==PositionNEER
        title('NEER')
    elseif xx==PositionRealGDP
        title('Real GDP')
    elseif xx==PositionPrices
        title('Prices')
    elseif xx==PositionRepo7Day
        title('Repo 7Day')
    else
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   Variance decompositions
FractionsOfVariance=zeros(Horizon+1,N,NNN);
xx=0;
hh=1;
while hh<=NumberOfDraws
    if A0IN(hh)==1
        xx=xx+1;
        FractionsOfVariance(:,:,xx)=GetFractionsOfVariance([zeros(N,1) BBBB(:,:,hh)],AA00(:,:,hh),N,LagOrder,Horizon,1);
    end
    hh=hh+1;
end
SortedFractionsOfVariance=sort(FractionsOfVariance,3);
HOR=(0:1:Horizon)';
%
Percentiles=fix(NNN*[0.5 0.16 0.84 0.05 0.95]');
%
PercentilesFractionsOfVariance=SortedFractionsOfVariance(:,:,Percentiles);
%
for xx=1:N
    Perc=squeeze(PercentilesFractionsOfVariance(:,xx,:));
    figure(1)
    subplot(2,N,N+xx)
    plot(HOR,Perc(:,1),'k',HOR,Perc(:,2:3),'r:',HOR,Perc(:,4:5),'r','LineWidth',2)
    axis([0 Horizon 0 1])      
    if xx==PositionRealInvestment
        title('Investment')
        ylabel('FEVs')
    elseif xx==PositionRealConsumption
        title('Consumption')
    elseif xx==PositionRealImport
        title('Import')
    elseif xx==PositionRealExport
        title('Export')
    elseif xx==PositionLogM2
        title('M2')
    elseif xx==PositionSpread
        title('Spread')
    elseif xx==PositionNEER
        title('NEER')
    elseif xx==PositionRealGDP
        title('Real GDP')
    elseif xx==PositionPrices
        title('Prices')
    elseif xx==PositionRepo7Day
        title('Repo 7Day')
    else
    end
end




