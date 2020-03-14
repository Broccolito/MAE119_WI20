% For solar Capacity
SurfTilt=10;
SurfAz=180;
tyadata=pvl_readtmy3('725905TYA.CSV');
TimeMatlab = tyadata.DateNumber;
Time = pvl_maketimestruct(TimeMatlab, ones(size(TimeMatlab))*tyadata.SiteTimeZone);
dayofyear = pvl_date2doy(Time.year, Time.month, Time.day);
DNI = tyadata.DNI; 
DHI = tyadata.DHI;
GHI = tyadata.GHI;
location = pvl_makelocationstruct(tyadata.SiteLatitude,tyadata.SiteLongitude,...
    tyadata.SiteElevation); %Altitude is optional
pressure= tyadata.Pressure*100; %Convert pressure from mbar to Pa
[SunAz, SunEl, ApparentSunEl, SolarTime] = pvl_ephemeris(Time, location);
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

Tamb=tyadata.DryBulb(1:8760);
windspeed=tyadata.Wspd(1:8760);
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

AC_power=pvl_snlinverter(Inverter,Vdc,Pdc);
AC_power(AC_power < 0) = 0;

mean_Jan = sum(reshape(AC_power(1:744),24,31),2)/31;
mean_Feb = sum(reshape(AC_power(745:1416),24,28),2)/28;
mean_Mar = sum(reshape(AC_power(1417:2160),24,31),2)/31;
mean_Apr = sum(reshape(AC_power(2161:2880),24,30),2)/30;
mean_May = sum(reshape(AC_power(2881:3624),24,31),2)/31;
mean_Jun = sum(reshape(AC_power(3625:4344),24,30),2)/30;
mean_Jul = sum(reshape(AC_power(4345:5088),24,31),2)/31;
mean_Aug = sum(reshape(AC_power(5089:5832),24,31),2)/31;
mean_Sep = sum(reshape(AC_power(5833:6552),24,30),2)/30;
mean_Oct = sum(reshape(AC_power(6553:7296),24,31),2)/31;
mean_Nov = sum(reshape(AC_power(7297:8016),24,30),2)/30;
mean_Dec = sum(reshape(AC_power(8017:8760),24,31),2)/31;

power_generation = [mean_Jan;mean_Feb;mean_Mar;mean_Apr;mean_May;mean_Jun;mean_Jul;mean_Aug;...
    mean_Sep;mean_Oct;mean_Nov;mean_Dec]*10^-6; %MW

fig = figure('units','inch','position',[5,5,6,5]);
hold on
plot(1:288, power_generation)
xlabel('Typical Hours of a year')
ylabel('Solar Power generated')
hold off
print(fig,'fig1.png','-dpng','-r800');