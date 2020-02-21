%% P1.2
DBfile = 'SandiaModuleDatabase_20120925.xlsx';
Module = pvl_sapmmoduledb(169, DBfile);

MS=28; % #module in series
MP=17; % # number of parallel strings 
Tamb=weather.DryBulb;
windspeed=weather.WindSpeed;
a=Module.a_wind;
b=Module.b_wind;
deltaT=Module.delT;



E0=1000;
SF=.98;


F1=max(0,polyval(Module.a,AM));
F2=max(0,polyval(Module.b,AOI));
Ee=F1.*((Eb.*F2+Module.fd * Ediff)/E0)*SF;
Ee(isnan(Ee))=0;


Tcell= pvl_sapmcelltemp(POA, 1000, a, b, windspeed, Tamb, deltaT);
Result = pvl_sapm(Module, Ee, Tcell);


VDC=MS*Result.Vmp;
IDC=MP*Result.Imp;
VDC(isnan(VDC))=0;
IDC(isnan(IDC))=0;
VDC(VDC<0)=0;

IDC(IDC<0)=0;
pwrDC=(VDC.*IDC);
yrlyDC_tot=2*sum(pwrDC) %includes two arrays
maxTcell=max(Tcell)

%% P1.3
load 'SandiaInverterDatabaseSAM2014.1.14.mat';
Inverter = SNLInverterDB(1393);
pwrAC= pvl_snlinverter(Inverter, VDC,pwrDC);
pwrAC(pwrAC<0)=0;
pwrAC_act=pwrAC; 
yrlyAC_tot=2*sum(pwrAC_act) %since 2 inverters

Inv_eff=pwrAC./pwrDC;
Inv_eff(isnan(Inv_eff))=0;
pwr_factor=pwrDC/Inverter.Pac0;
[maxDC ind_maxDC]=max(pwrDC);
[minDC ind_minDC]=min(pwrDC);

month_DCmax=Time.month(ind_maxDC);
day_DCmax=Time.day(ind_maxDC);

month_DCmin=Time.month(ind_minDC);
day_DCmin=Time.day(ind_minDC);
% 
figure
hold on
tfilterMAXDC= and(Time.month == month_DCmax,Time.day == day_DCmax);
tfilterMINDC= and(Time.month == month_DCmin,Time.day == day_DCmin);
plot(Time.hour(tfilterMAXDC),pwrDC(tfilterMAXDC),'-s')
plot(Time.hour(tfilterMAXDC),pwrAC_act(tfilterMAXDC),'-s')
title('DC VS AC power during largest daily DC total');
xlabel('Hour of the Day');
ylabel('Power (W)');
legend('DC Power', 'AC Power');

figure
hold on
tfilterMINDC= and(Time.month == month_DCmin,Time.day == day_DCmin);
plot(Time.hour(tfilterMINDC),pwrDC(tfilterMINDC),'-s')
plot(Time.hour(tfilterMINDC),pwrAC_act(tfilterMINDC),'-s')
title('DC VS AC power during smallest daily DC total');
xlabel('Hour of the Day');
ylabel('Power (W)');
legend('DC Power', 'AC Power');

figure
scatter(Inv_eff,pwr_factor)
title('Inverter Efficiency vs. Inverter Power Factor')
xlabel('Inverter Efficiency')
ylabel('Inverter Power Factor')