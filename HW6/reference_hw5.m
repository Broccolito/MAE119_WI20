% MAE 119 HW 5
% Claire Stone A13947716
% Trish Tran A13710432

clear;clc;
addpath(genpath('pvlib'));

%Load demand
dem = readtable('Konocti_Load.csv');
Demand = table2array(dem(:,3)); % MW

%Panel and inverter configuration
%Panel: Canadian solar CS5P-220M
ModuleParameters = pvl_sapmmoduledb(123,'SandiaModuleDatabase_20120925.xlsx');
%Inverter: Green Power Technologies PV900WD
load('SandiaInverterDatabaseSAM2014.1.14.mat')
Inverter = SNLInverterDB(441);
clear InverterNames SNLInverterDB

%Set up array configuration
Array.Tilt = 38;
Array.Azimuth = 180; %south
ArrayMs = 100; %Number of modules in series
ArrayMp = 50; %Number of parallel strings

% Find solar generation over the year
%Set location and time - using 2019
Location = pvl_makelocationstruct(39,-123,192); %Ukiah,CA lat, lon and alt
times = datetime(2019,1,1,0,0,0):hours(1):datetime(2019,12,31,23,0,0); %2019 with 1h resolution
Time = pvl_maketimestruct(datenum(times),-8); % Feed in times and time zone

n = 1:365;
doy = repelem(n,24);

%Load TMY3 data for Ukiah
TMYData = pvl_readtmy3('725905TYA.CSV'); %Ukiah, CA TMY3
%Get DNI, DIF, and GHI
DNI = TMYData.DNI;
DIF = TMYData.DHI;
GHI = TMYData.GHI;

%Get realistic angles
PresPa = TMYData.Pressure*100; %Pressure must be in Pascals
DryBulb = TMYData.DryBulb;
[SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location,PresPa,DryBulb);
AppZ = 90-AppSunEl;

%Get angle of incidence
AOI = pvl_getaoi(Array.Tilt, Array.Azimuth, AppZ, SunAz);

% Calculate POA
%Direct component
Eb = 0*AOI; %Initialize variable
Eb(AOI<90) = DNI(AOI<90).*cosd(AOI(AOI<90));

%Diffuse component
HExtra = pvl_extraradiation(doy);
AM = pvl_relativeairmass(AppZ); %air mass
%Get sky diffuse from Perez model
EdiffSky = pvl_perez(Array.Tilt, Array.Azimuth, DIF, DNI, HExtra, AppZ, SunAz, AM);
Albedo = 0.2;
%Get ground reflected diffuse
GHI(isnan(GHI))= 0;
EdiffGround = pvl_grounddiffuse(Array.Tilt, GHI, Albedo);

% Arrange POA
E = Eb + EdiffSky + EdiffGround; %Total incident irradiance (W/m^2)
Ediff = EdiffSky + EdiffGround; %Total diffuse incident irradiance (W/m^2)

%Get PV cell operating temperatures
SF = 0.98; %2 percent soiling level
E0 = 1000; %Reference irradiance
%Get cell temp
celltemp = pvl_sapmcelltemp(E, E0, ModuleParameters.a_wind,...
    ModuleParameters.b_wind, TMYData.Wspd, DryBulb, ModuleParameters.delT);

% Correct POA (spectral losses)
F1 = max(0,polyval(ModuleParameters.a,AM)); %Spectral loss function
F2 = max(0,polyval(ModuleParameters.b,AOI)); %Angle of incidence loss function
Ee = F1.*((Eb.*F2+ModuleParameters.fd.*Ediff)/E0)*SF; %Effective irradiance
Ee(isnan(Ee)) = 0; %Set any NaNs to zero

% Use Sandia Array Performance Model to get DC power
mSAPMResults = pvl_sapm(ModuleParameters, Ee, celltemp);
mSAPMResults.Vmp(mSAPMResults.Vmp<0)=0; %set any spurious negatives to zero
mSAPMResults.Imp(mSAPMResults.Imp<0)=0; %set any spurious negatives to zero
Vdc = ArrayMs*mSAPMResults.Vmp; %Voltage summed by panels in series
Idc = ArrayMp*mSAPMResults.Imp; %Current summed by panels in parallel
Pdc_1array = Vdc .*Idc; %Total 1 array power vector

%Get AC power using the inverter specs
Pac = pvl_snlinverter(Inverter,Vdc,Pdc_1array);
Pac(Pac<0)=0; %set any spurious negatives to zero

Pac_Jan = reshape(Pac(1:744),24,31);
Pac_Feb = reshape(Pac(745:1416),24,28);
Pac_Mar = reshape(Pac(1417:2160),24,31);
Pac_Apr = reshape(Pac(2161:2880),24,30);
Pac_May = reshape(Pac(2881:3624),24,31);
Pac_Jun = reshape(Pac(3625:4344),24,30);
Pac_Jul = reshape(Pac(4345:5088),24,31);
Pac_Aug = reshape(Pac(5089:5832),24,31);
Pac_Sep = reshape(Pac(5833:6552),24,30);
Pac_Oct = reshape(Pac(6553:7296),24,31);
Pac_Nov = reshape(Pac(7297:8016),24,30);
Pac_Dec = reshape(Pac(8017:8760),24,31);

AvgDayJan = sum(Pac_Jan,2)/31;
AvgDayFeb = sum(Pac_Feb,2)/28;
AvgDayMar = sum(Pac_Mar,2)/31;
AvgDayApr = sum(Pac_Apr,2)/30;
AvgDayMay = sum(Pac_May,2)/31;
AvgDayJun = sum(Pac_Jun,2)/30;
AvgDayJul = sum(Pac_Jul,2)/31;
AvgDayAug = sum(Pac_Aug,2)/31;
AvgDaySep = sum(Pac_Sep,2)/30;
AvgDayOct = sum(Pac_Oct,2)/31;
AvgDayNov = sum(Pac_Nov,2)/30;
AvgDayDec = sum(Pac_Dec,2)/31;

SolarOutput = [AvgDayJan;AvgDayFeb;AvgDayMar;AvgDayApr;AvgDayMay;AvgDayJun;AvgDayJul;AvgDayAug;...
    AvgDaySep;AvgDayOct;AvgDayNov;AvgDayDec]*10^-6; %MW

%initialize LCC matrix
LCC = zeros(200,80);

for k = 10:80 %solar capacity, MW
    for j = 10:200 %battery capacity, MWh
        
        SOL = SolarOutput*k; %solar generation * number of arrays, MW

        Battery = zeros(288,1); %record charge of battery each hour
        Battery(1) = j;
        PowerOut = zeros(288,1); %record number of unserved MW (will sum later)

        for i = 1:288   
               if SOL(i) < Demand(i)
                    Battery(i+1) = Battery(i) - (Demand(i)-SOL(i)); %charge battery
                    if Battery(i+1) < 0
                        PowerOut(i) = abs(Battery(i+1)); %record # of unserved MWh
                        Battery(i+1) = 0; %reset battery to zero
                    end
                            
                elseif SOL(i) >= Demand(i)
                    Battery(i+1) = Battery(i) + (SOL(i) - Demand(i));
                    if Battery(i+1) > j
                       Battery(i+1) = j;
                    end   
               end
        end

        % LCC analysis - costs are positive, gains are negative
        
        %Initial vectors
        PO = zeros(21,1); %power outtage cost
        Invest = zeros(21,1); %Investment costs for solar and battery
        OM = zeros(21,1); %O&M costs for solar
        Incent = zeros(21,1); %incentives for solar

        SCap_kw = k*1000; %solar capacity in kW
        BCap_kw = j*1000; %battery capacity in kW

        PO(1) = sum(PowerOut)*1000*14; % $14/unserved kWh
        Invest(1) = 1800*SCap_kw + 600*BCap_kw; %solar and battery investment cost in Y0
        Invest(11) = 600*BCap_kw; %new battery investment cost after 10 yr
        OM(2) = 9*SCap_kw; % O&M for solar
        Incent(1) = -1800*SCap_kw*0.3;
        

        %inflation for O&M
        for m = 3:21
            OM(m) = OM(2)*(1.05)^(m-2);
        end
        
        %Inflation for power outage costs
        for p = 2:21
            PO(p) = PO(1)*(1.05)^(p-2); 
        end
        
        %Sum costs by year
        sumcost = Invest + PO + OM + Incent;
        
        %discount each year
        discount = zeros(21,1);
        for n = 1:21
            discount(n) = sumcost(n)/(1.06)^(n-1);
        end
        
        %Sum the discounted totals to get LCC. store in matrix
        LCC(j,k) = sum(discount);
    end
end

% %remove zeros from LCC matrix
LCC = LCC(10:end,10:end);
X = 10:80;
Y = 10:200;

figure(1)
surf(X,Y,LCC)

figure(2)
contourf(X,Y,LCC,'LineStyle','none')
colorbar

lowest = min(min(LCC));
[b,s] = find(LCC == lowest);

bestbattery = Y(b);
bestsolar = X(s);

fprintf('The best combination is a %d MW solar system and a %d MWh battery.\n',bestsolar,bestbattery)
fprintf('This has a 20 yr LCC of $%d.\n\n',lowest)

fprintf('In comparison, the ad-hoc sizing of HW 4 was a 34 MW solar system and an 85 MWh battery, with an LCC of $127,508,920.\n')
fprintf('The ad-hoc sizing was smaller because we only calculated it for 10 days usage in October.\n')

