clear;
close all;
warning off
addpath(genpath('PVLib'))

power_load = readtable('Konocti_Load.csv');

% For solar Capacity
SurfTilt=10;
SurfAz=180;
TMYData=pvl_readtmy3('725905TYA.CSV'); %Ukiah, CA site
TimeMatlab = TMYData.DateNumber;
Time = pvl_maketimestruct(TimeMatlab, ones(size(TimeMatlab))*TMYData.SiteTimeZone);
dayofyear = pvl_date2doy(Time.year, Time.month, Time.day);
DNI = TMYData.DNI; % Read in for comparison with results
DHI = TMYData.DHI; % Read in for comparison with results
GHI = TMYData.GHI;
Location = pvl_makelocationstruct(TMYData.SiteLatitude,TMYData.SiteLongitude,...
    TMYData.SiteElevation); %Altitude is optional
pressure= TMYData.Pressure*100; %Convert pressure from mbar to Pa
[SunAz, SunEl, ApparentSunEl, SolarTime] = pvl_ephemeris(Time, Location);
SunZen=90-ApparentSunEl;
AM= pvl_relativeairmass(SunZen);
AMa=pvl_absoluteairmass(AM,pressure);

HExtra = pvl_extraradiation(dayofyear);
Ediffsky = pvl_perez(SurfTilt, SurfAz, DHI, DNI, HExtra, SunZen, SunAz, AMa);

ro_g = 0.2;
AOI = pvl_getaoi(SurfTilt, SurfAz, SunZen, SunAz);
Eb=0*AOI;
Eb(AOI<90)=DNI(AOI<90).*cosd(AOI(AOI<90));
GHI(isnan(GHI))=0;
Ediffground=pvl_grounddiffuse(SurfTilt, GHI, ro_g);

%POA
POA=Eb + Ediffsky + Ediffground;
Ediff=Ediffsky+Ediffground;

%Solar panel
DBfile = 'SandiaModuleDatabase_20120925.xlsx';
Module = pvl_sapmmoduledb(124, DBfile);

Tamb=TMYData.DryBulb(1:8760);
windspeed=TMYData.Wspd(1:8760);
a=Module.a_wind;
b=Module.b_wind;
deltaT=Module.delT;

Ee=POA*.98/1000;
Tcell = pvl_sapmcelltemp(Ee, 1000, a, b, windspeed, Tamb, deltaT);
Result = pvl_sapm(Module, Ee, Tcell);

MS=100; % #module in series
MP=50; % # number of parallel strings

Vdc=MS*Result.Vmp;
Vdc(Vdc<0)=0;
Idc=MP*Result.Imp;
Pdc=(Vdc.*Idc);
load 'SandiaInverterDatabaseSAM2014.1.14.mat';
Inverter = SNLInverterDB(441);

Pac=pvl_snlinverter(Inverter,Vdc,Pdc);
Pac(Pac < 0) = 0;
Oct_solar_generation = sum(Pac);
max_solarPwr=max(Pac);

PacJan = Pac(337:360);
PacFeb = Pac(1081:1104);
PacMar = Pac(1753:1776);
PacApr = Pac(2497:2520);
PacMay = Pac(3217:3240);
PacJun = Pac(3961:3984);
PacJul = Pac(4681:4704);
PacAug = Pac(5425:5448);
PacSep = Pac(6169:6192);
PacOct = Pac(6889:6912);
PacNov = Pac(7633:7656);
PacDec = Pac(8353:8376);

fig = figure('units','inch','position',[5,5,6,5]);
hold on
plot(1:24, PacJan)
plot(1:24, PacFeb)
plot(1:24, PacMar)
plot(1:24, PacApr)
plot(1:24, PacMay)
plot(1:24, PacJun)
plot(1:24, PacJul)
plot(1:24, PacAug)
plot(1:24, PacSep)
plot(1:24, PacOct)
plot(1:24, PacNov)
plot(1:24, PacDec)
legend('Jan', 'Feb',	'Mar',	'Apr',	'May',	'Jun',	'Jul',...
    'Aug',	'Sep',	'Oct',	'Nov',	'Dec')
xlabel('Hour of a typical day')
ylabel('Solar Power generated')
hold off
print(fig,'fig2.png','-dpng','-r800');

battery_threshold = 5e6;
generation = [PacJan;PacFeb;PacMar;PacApr;PacMay;PacJun;PacJul;PacAug;PacSep;PacOct;PacNov;PacDec];
load_generation = [power_load, table(generation)];
load_generation.Konocti_Load = load_generation.Konocti_Load * 1e5;
supply = load_generation.generation - load_generation.Konocti_Load;
residual = zeros(288,1);
residual(1) = battery_threshold; % Start with full battery
for i = 2:288
    if residual <= battery_threshold
        residual(i) = residual(i-1) + supply(i);
        if residual(i) < 0
            residual(i) = 0;
        end
    else
        residual(i) = residual(i-1);
    end
end

fig = figure('units','inch','position',[5,5,6,5]);
hold on
plot(residual)
xlabel('Hour with average power generation profile')
ylabel('Power left in the battery (w)')
legend('When power left is zero, outage occurs')
hold off
print(fig,'fig2.png','-dpng','-r800');

outage_day_count = sum((residual == 0));
