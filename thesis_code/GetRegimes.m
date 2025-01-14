function Regimes=GetRegimes(Time,Dates)
for hh=1:length(Dates)
    Regimes(:,hh)=(Time>=Dates(hh))*(-1)^hh;
end
if length(Dates)>1
    Regimes=(sum(Regimes')'*2+1)*1000000;
else
    Regimes=(Regimes*2+1)*1000000;
end
