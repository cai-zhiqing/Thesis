%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
path(path,'C:\EmpiricalMacro')
options = optimset;
options = optimset(options,'Display','off');
warning('off','all')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumberOfDraws=2000;
LagOrder=6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       Loading time series:
X=xlsread('C:\EmpiricalMacro\MonetaryPolicyHousePrices.xls','USMonthly','C157:V738');
[T,N]=size(X);
Time=X(:,1); % Column 1: a numerical string representing calendar time. E.g., 2010Q1 is 2010.25, 2010 Q2 is 2010.5, etc.
LogHousePrice=log(X(:,2));
LogRent=log(X(:,3));
LogCPI=log(X(:,4));
UnemploymentRate=X(:,5)/100;
FEDFUNDS=X(:,7)/100;
LogHousingStarts=log(X(:,6));
LogIndustrialProduction=log(X(:,19));
LogNonBorrowedReserves=log(X(:,9));
LogTotalReserves=log(X(:,10));
LogMarkWatsonRealGDP=log(X(:,11));
LogPCEDeflator=log(X(:,13));
LogCommoditiesPriceIndex=log(X(:,14));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=[LogNonBorrowedReserves LogTotalReserves LogCommoditiesPriceIndex UnemploymentRate LogMarkWatsonRealGDP LogPCEDeflator FEDFUNDS];
[T,N]=size(X);
Stationary='N';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dates=[1971+8/12 1979+8/12 2007+8/12 2008+10/12]';
Regimes=GetRegimes(Time,Dates);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Trend=0;
Intercept='Y';
Horizon=5*12;
T=length(Time);
NN=NumberOfDraws;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randn('state',1234)
rand('state',1234)
%
[BBBB,SSSS,UUUU]=BayesianVARp(X,LagOrder,Intercept,NumberOfDraws,Stationary);
if Intercept=='Y'
    MUMU=BBBB(:,1,:);
    BBBB=BBBB(:,2:size(BBBB,2),:);
else
    MUMU=zeros(N,1,NN);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear DFILE varname
DFILE=['C:\EmpiricalMacro\BayesianVARReducedFormEstimates_ACRR'];
varname(1,:)=['BBBB'];
varname(2,:)=['MUMU'];
varname(3,:)=['SSSS'];
varname(4,:)=['UUUU'];
%
save(DFILE,varname(1,:),varname(2,:),varname(3,:),varname(4,:))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              What follows is based on combining inertial and sign restrictions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear DFILE varname
DFILE=['C:\EmpiricalMacro\A0Matrix_ACRR'];
varname(1,:)=['AA00'];
varname(2,:)=['A0IN'];
%
AA00=zeros(N,N,NumberOfDraws);
A0IN=zeros(NumberOfDraws,1);
%
randn('state',1234)
rand('state',1234)
%
hh=1;
while hh<=NumberOfDraws
    NumberOfDraws-hh+1
    var=SSSS(:,:,hh);
    [A0,IND]=GetA0AriasCaldaraRubioRamirez(var);
    if IND==1
        AA00(:,:,hh)=A0;
        A0IN(hh)=IND;
    end
    hh=hh+1;
end
save(DFILE,varname(1,:),varname(2,:))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NNN=sum(sum(A0IN)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         II: Impulse-response functions
IRFs=zeros(Horizon+1,N,NNN);
%
xx=0;
hh=1;
while hh<=NumberOfDraws
    % NumberOfDraws-hh+1
    bb=BBBB(:,:,hh);
    %
    if A0IN(hh)==1
        xx=xx+1;
        A0=AA00(:,:,hh);
        IRFs(:,:,xx)=GetIRFs(bb,A0,N,LagOrder,Horizon,1);
    end
    hh=hh+1;
end
%
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
    if xx==1
        title('Non borrowed reserves')
        ylabel('IRFs')
    elseif xx==2
        title('Total reserves')
    elseif xx==3
        title('Commodities price index')
    elseif xx==4
        title('Unemployment rate')
    elseif xx==5
        title('Real GDP')
    elseif xx==6
        title('PCE deflator')
    elseif xx==7
        title('Federal Funds rate')
    else
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              III: Variance decompositions
clear DFILE varname
DFILE=['C:\EmpiricalMacro\FractionsOfVariance_ACRR'];
varname(1,:)=['FractionsOfVariance'];
%
FractionsOfVariance=zeros(Horizon+1,N,NNN);
%
xx=0;
hh=1;
while hh<=NumberOfDraws
    NumberOfDraws-hh+1
    mu=MUMU(:,:,hh);
    bb=BBBB(:,:,hh);
    if A0IN(hh)==1
        xx=xx+1;
        A0=AA00(:,:,hh);
        FractionsOfVariance(:,:,xx)=GetFractionsOfVariance([mu bb],A0,N,LagOrder,Horizon,1);
    end
    hh=hh+1;
end
save(DFILE,varname(1,:))
%
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
    if xx==1
        title('Non borrowed reserves')
        ylabel('Fractions of FEV')
    elseif xx==2
        title('Total reserves')
    elseif xx==3
        title('Commodities price index')
    elseif xx==4
        title('Unemployment rate')
    elseif xx==5
        title('Real GDP')
    elseif xx==6
        title('PCE deflator')
    elseif xx==7
        title('Federal Funds rate')
    else
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                           IV: Counterfactuals killing off the shocks
CounterfactualX=zeros(size(X,1),size(X,2),NNN);
CounterfactualMinusActualX=zeros(size(X,1),size(X,2),NNN);
%
xx=0;
hh=1;
while hh<=NumberOfDraws
    NumberOfDraws-hh+1
    mu=MUMU(:,:,hh);
    bb=BBBB(:,:,hh);
    uu=UUUU(:,:,hh);
    %
    if A0IN(hh)==1
        xx=xx+1;
        A0=AA00(:,:,hh);
        shocks=A0\uu;
        shocks(1,:)=zeros(size(shocks(1,:)));
        counterfactualx=X';
        for tt=LagOrder+1:size(X,1)
            counterfactualx(:,tt)=mu+bb*vec(fliplr(counterfactualx(:,tt-LagOrder:tt-1)))+A0*shocks(:,tt-LagOrder);
        end
        counterfactualx=counterfactualx';
        %
        CounterfactualX(:,:,xx)=counterfactualx;
        CounterfactualMinusActualX(:,:,xx)=counterfactualx-X;
    end
    hh=hh+1;
end
%
for xx=1:N
    PercentilesCounterfactual=ExtractPercentiles(sort(squeeze(CounterfactualX(:,xx,:))'),[0.5 0.16 0.84 0.05 0.95]')';
    figure(2)
    subplot(2,N,xx)
    plot(Time,X(:,xx),'b',Time,PercentilesCounterfactual(:,1),'k',Time,PercentilesCounterfactual(:,2:3),'r')
    xlim([Time(LagOrder+1) Time(T)])
    if xx==1
        title('Non borrowed reserves')
    elseif xx==2
        title('Total reserves')
    elseif xx==3
        title('Commodities price index')
    elseif xx==4
        title('Unemployment rate')
    elseif xx==5
        title('Real GDP')
    elseif xx==6
        title('PCE deflator')
    elseif xx==7
        title('Federal Funds rate')
    else
    end
    PercentilesCounterfactualMinusActualX=ExtractPercentiles(sort(squeeze(CounterfactualMinusActualX(:,xx,:))'),[0.5 0.16 0.84 0.05 0.95]')';
    figure(2)
    subplot(2,N,N+xx)
    plot(Time,zeros(size(Time)),'b:',Time,PercentilesCounterfactualMinusActualX(:,1),'k',Time,PercentilesCounterfactualMinusActualX(:,2:3),'r')
    xlim([Time(LagOrder+1) Time(T)])
end
%
disp('---------------------------------------------------------------------------------')
disp('So you see that killing off monetary policy shocks not much happens: these shocks')
disp('these shocks do not play much of a role for anything ...')
disp(' ')
disp('Now lets see a counterfactual in which we make the structural monetary policy')
disp('rule more aggressive against inflation ...')
disp('---------------------------------------------------------------------------------')
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %                                          V: A counterfactual in which we make the structural monetary
% % % %                                               policy rule more aggressive against the PCE deflator
% % % K=1.05; % This is the factor by which we increase the aggressiveness against the PCE deflator
% % % %
% % % CounterfactualX=zeros(size(X,1),size(X,2),NNN);
% % % CounterfactualMinusActualX=zeros(size(X,1),size(X,2),NNN);
% % % %
% % % INDEX=zeros(NNN,1);
% % % %
% % % xx=0;
% % % hh=1;
% % % while hh<=NumberOfDraws
% % %     NumberOfDraws-hh+1
% % %     mu=MUMU(:,:,hh);
% % %     bb=BBBB(:,:,hh);
% % %     uu=UUUU(:,:,hh);
% % %     %
% % %     Index=find(A0IN(:,hh));
% % %     zz=1;
% % %     while zz<=length(Index)
% % %         xx=xx+1;
% % %         A0=AA00(:,:,Index(zz),hh);
% % %         shocks=A0\uu;
% % %         %
% % %         InvA0=MyInverse(A0);
% % %         StructuralMU=InvA0*mu;
% % %         StructuralB=InvA0*bb;
% % %         StructuralBOld=StructuralB;
% % %         InvA0(1,N-1)=K*InvA0(1,N-1);
% % %         for ll=1:LagOrder
% % %             StructuralB(1,(ll-1)*N+(N-1))=K*StructuralB(1,(ll-1)*N+(N-1));
% % %         end
% % %         A0=MyInverse(InvA0);
% % %         MU=A0*StructuralMU;
% % %         B=A0*StructuralB;
% % %         %
% % %         R=varroots(LagOrder,N,[MU B]);
% % %         if max(abs(R))<1
% % %             INDEX(xx)=1;
% % %             counterfactualx=X';
% % %             for tt=LagOrder+1:size(X,1)
% % %                 counterfactualx(:,tt)=MU+B*vec(fliplr(counterfactualx(:,tt-LagOrder:tt-1)))+A0*shocks(:,tt-LagOrder);
% % %             end
% % %             counterfactualx=counterfactualx';
% % %             %
% % %             CounterfactualX(:,:,xx)=counterfactualx;
% % %             CounterfactualMinusActualX(:,:,xx)=counterfactualx-X;
% % %         end
% % %         zz=zz+1;
% % %     end
% % %     hh=hh+1;
% % % end
% % % %
% % % Indices=find(INDEX);
% % % CounterfactualX=CounterfactualX(:,:,Indices);
% % % CounterfactualMinusActualX=CounterfactualMinusActualX(:,:,Indices);
% % % %
% % % for xx=1:N
% % %     PercentilesCounterfactual=ExtractPercentiles(sort(squeeze(CounterfactualX(:,xx,:))'),[0.5 0.16 0.84 0.05 0.95]')';
% % %     figure(3)
% % %     subplot(2,N,xx)
% % %     plot(Time,X(:,xx),'b',Time,PercentilesCounterfactual(:,1),'k',Time,PercentilesCounterfactual(:,2:3),'r')
% % %     xlim([Time(LagOrder+1) Time(T)])
% % %     if xx==1
% % %         title('Non borrowed reserves')
% % %     elseif xx==2
% % %         title('Total reserves')
% % %     elseif xx==3
% % %         title('Commodities price index')
% % %     elseif xx==4
% % %         title('Unemployment rate')
% % %     elseif xx==5
% % %         title('Real GDP')
% % %     elseif xx==6
% % %         title('PCE deflator')
% % %     elseif xx==7
% % %         title('Federal Funds rate')
% % %     else
% % %     end
% % %     PercentilesCounterfactualMinusActualX=ExtractPercentiles(sort(squeeze(CounterfactualMinusActualX(:,xx,:))'),[0.5 0.16 0.84 0.05 0.95]')';
% % %     figure(3)
% % %     subplot(2,N,N+xx)
% % %     plot(Time,zeros(size(Time)),'b:',Time,PercentilesCounterfactualMinusActualX(:,1),'k',Time,PercentilesCounterfactualMinusActualX(:,2:3),'r')
% % %     xlim([Time(LagOrder+1) Time(T)])
% % % end
% % % %
% % % disp('-----------------------------------------------------------------------------------------')
% % % disp('So you see that if you increase the aggressiveness of the structural monetary policy rule')
% % % disp('against inflation by a sufficient amount, you get a lower inflation path ...')
% % % disp('     ')
% % % disp('The interesting thing to notice is that this does not seem to have a big impact on either')
% % % disp('the unemployment rate or real GDP growth ...')
% % % disp('-------------------------------------------------------------------------------')
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
