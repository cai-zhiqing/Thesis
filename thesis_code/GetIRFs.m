function IRFs=GetIRFs(B,A0,N,LagOrder,Horizon,K)
%
YY=zeros(N,LagOrder+1+Horizon);
Shocks=zeros(N,1);
Shocks(K)=1;
YY(:,LagOrder+1)=A0*Shocks;
for tt=LagOrder+2:LagOrder+1+Horizon
    YY(:,tt)=B*vec(fliplr(YY(:,tt-LagOrder:tt-1)));
end
YY=YY';
YY=YY(LagOrder+1:LagOrder+1+Horizon,:);
IRFs=YY;








