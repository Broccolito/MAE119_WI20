%% Homework 4
% Reuben Lumaban
% Eni Ikuku

%% P1.1
clear all;
close all;
clc;
addpath(genpath('pvlib'));
%for a typical electrical demand on day in Oct
hour = linspace(0,23,24);
demand=[5.1 4.9 5.1 5.3 5.6 6.7 8.0 8.6 8.3 7.6 7.2 6.7 6.2 5.9 5.9 5.9 6.5 7.2 7.9 8.1 8.1 7.4 6.4 5.6]*10^6; %W


peak_d=max(demand); % peak demand

% Nominal Capacity (Monthly)-31 days in Oct
diesel_nomcap=peak_d % Diesel ; eql to peak load on a typical day

btry=[sum(demand(20:end)) sum(demand(1:8))]; %Battery ; nighttime energy use 
battery_nomcap=max(btry)

% For solar Capacity
SurfTilt=10;
SurfAz=180;
TMYData=pvl_readtmy3('725905TYA.csv'); %Ukiah, CA site
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

ro_g=.2;
AOI = pvl_getaoi(SurfTilt, SurfAz, SunZen, SunAz);
Eb=0*AOI;
Eb(AOI<90)=DNI(AOI<90).*cosd(AOI(AOI<90));
GHI(isnan(GHI))=0;
Ediffground=pvl_grounddiffuse(SurfTilt, GHI, ro_g);

% For DOY:
DOYjan = 1:31;
DOYfeb = 32:59;
DOYmar = 60:90;
DOYapr = 91:120;
DOYmay = 121:151;
DOYjun = 152:181;
DOYjul = 182:212;
DOYaug = 213:243;
DOYsept = 244:273;
DOYoct = 274:304;
DOYnov = 305:334;
DOYdec = 335:365;

jan = DOYjan*24;
feb = DOYfeb*24;
mar = DOYmar*24;
apr = DOYapr*24;
may = DOYmay*24;
jun = DOYjun*24;
jul = DOYjul*24;
aug = DOYaug*24;
sept = DOYsept*24;
oct = DOYoct*24;
nov = DOYnov*24;
dec = DOYdec*24;

%POA
POA=Eb + Ediffsky + Ediffground;
Ediff=Ediffsky+Ediffground;

POA=POA(oct); %actual generation in Oct

%monthly demand
Oct_demand = sum(demand)*31; %W


%Solar panel
DBfile = 'SandiaModuleDatabase_20120925.xlsx';
Module = pvl_sapmmoduledb(124, DBfile);


Tamb=TMYData.DryBulb(oct);
windspeed=TMYData.Wspd(oct);
a=Module.a_wind;
b=Module.b_wind;
deltaT=Module.delT;

Ee=POA*.98/1000;
Tcell = pvl_sapmcelltemp(Ee, 1000, a, b, windspeed, Tamb, deltaT);
Result = pvl_sapm(Module, Ee, Tcell);

%change MS and MP to get slr to equal Oct_damand, 5000 panels
MS=100; % #module in series
MP=50; % # number of parallel strings 

Vdc=MS*Result.Vmp;
Vdc(Vdc<0)=0;
Idc=MP*Result.Imp;
Pdc=(Vdc.*Idc);
load 'SandiaInverterDatabaseSAM2014.1.14.mat';
Inverter = SNLInverterDB(441);

Pac=pvl_snlinverter(Inverter,Vdc,Pdc);
Pac(Pac<0)=0;
Oct_solar_generation=sum(Pac)
max_solarPwr=max(Pac);

Solar_array_nomcap=10^6
Solar_arrays_req=ceil((Oct_demand/Oct_solar_generation))

%%
%Problem 2
 % System Operation for 20 days
 
 %Solar generation October 15, n=15 in october vectors
Oct15I1 = 24*(15-1)+1;
Oct15I2 = 24*15;

Oct15ACGen = Pac(Oct15I1:Oct15I2)*30; 
% the hourly AC power for 30 arrays on october 15

timeOct15 = datetime(2019,10,15,0,0,0):hours(1):datetime(2019,10,15,23,0,0);
%October 15 2019 every 1 hour

%Plot demand and solar array output
figure(1)
plot(timeOct15,Oct15ACGen,timeOct15,demand)

matrix = [demand; Oct15ACGen'];

solar = min(matrix); 
%takes the min of each column which is the values of the solar curve except where it goes above demand

diesel = demand - solar;

DieselTwenty = sum(diesel)*20/10^3; %energy needed from diesel for 20 days in Solar + Diesel System
Dieselcap = max(diesel)/10^3; %kW, system size


%Last Part: Calculating how much extra solar energy we have.
%Find avg daily energy
AvgEnergy = sum(Pac)/31; %Wh/day

%Total unused energy per year
YrExtraE = AvgEnergy*(365-20)/10^3; %kWh
