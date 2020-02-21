clc; clear; close all;
addpath('PVLib\')

%% Question 1

GHI_struct = load('measuredGHI2016.mat');
pressure_struct = load('measuredMeteorological2016.mat');
GHI = GHI_struct.GHI;
time_GHI = GHI_struct.time_GHI;

Location = pvl_makelocationstruct(33,-117,137); % San Diego
Time = pvl_maketimestruct(datenum(time_GHI), -8);

[SunAz, SunEl, ApparentSunEl, SolarTime] = pvl_ephemeris(Time, Location);
Z = 90 - ApparentSunEl;

doy = pvl_date2doy(Time.year,Time.month,Time.day);

[DNI, DHI, Kt] = pvl_erbs(GHI', Z, doy);

SurfTilt = 10;
SurfAz = 180;
HExtra = pvl_extraradiation(pvl_date2doy(Time.year,Time.month,Time.day));
pressure = pressure_struct.Pressure;
AMa = pvl_absoluteairmass(pvl_relativeairmass(Z), pressure);
[SkyDiffuse,SkyDiffuse_Iso,SkyDiffuse_Cir,SkyDiffuse_Hor] = pvl_perez(SurfTilt, SurfAz, DHI, DNI, HExtra, Z, SunAz, AMa);

AOI = pvl_getaoi(SurfTilt, SurfAz, Z, SunAz);
Eb = 0*AOI;
Eb(AOI<90) = DNI(AOI<90).*cosd(AOI(AOI<90));
Albedo=0.2;
GHI(isnan(GHI)) = 0;
EdiffGround = pvl_grounddiffuse(SurfTilt, GHI, Albedo);
POA = Eb + SkyDiffuse + EdiffGround;
Ediff = SkyDiffuse + EdiffGround;

daily_POA = zeros(366,1);
for i = 1:366
    daily_POA(i) = sum(POA(doy == i));
end

daily_GHI = zeros(366,1);
for i = 1:366
    daily_GHI(i) = sum(GHI(doy == i));
end

day_min = find(daily_POA == min(daily_POA));
day_max = find(daily_POA == max(daily_POA));

fig = figure('units','inch','position',[5,5,5,5]);
hold on
plot((1:366), daily_POA, 'k-')
plot(day_min, daily_POA(day_min), 'ro')
plot(day_max, daily_POA(day_max), 'ro')
xlabel('Day of the Year')
ylabel('Irradiance (W·m^{-2})')
xlim([0,366])
legend('POA', 'Maxmium and Minimum', 'Location','NorthWest')
hold off
print(fig,'Q1.1.png','-dpng','-r800');

max_GHI = GHI(doy == day_max)';
max_POA = POA(doy == day_max);

fig = figure('units','inch','position',[5,5,5,5]);
hold on
plot(max_POA, 'bo-', 'MarkerFaceColor', 'b')
plot(max_GHI, 'go-', 'MarkerFaceColor', 'g')
xlabel('Hour of the Day')
ylabel('Irradiance (W·m^{-2})')
xlim([0,24])
legend('Max POA', 'Max GHI')
title('GHI and POA Plot at Day 296 (Maximum POA)')
hold off
print(fig,'Q1.2.png','-dpng','-r800');

min_GHI = GHI(doy == day_min)';
min_POA = POA(doy == day_min);

fig = figure('units','inch','position',[5,5,5,5]);
hold on
plot(min_POA, 'bo-', 'MarkerFaceColor', 'b')
plot(min_GHI, 'go-', 'MarkerFaceColor', 'g')
xlabel('Hour of the Day')
xlim([0,24])
legend('Min POA', 'Min GHI')
title('GHI and POA Plot at Day 351 (Minimum POA)')
hold off
print(fig,'Q1.3.png','-dpng','-r800');

%% Question 2
ModuleParameters = pvl_sapmmoduledb(169, 'SandiaModuleDatabase_20120925.xlsx');
E = Eb + SkyDiffuse + EdiffGround;
SF=0.98; % use 2% soiling levels
E0 = 1000; %Reference irradiance
celltemp = pvl_sapmcelltemp(E, E0, ModuleParameters.a_wind, ...
    ModuleParameters.b_wind, pressure_struct.WindSpeed, pressure_struct.DryBulb, ModuleParameters.delT);
F1 = max(0,polyval(ModuleParameters.a,AMa)); %Spectral loss function
F2 = max(0,polyval(ModuleParameters.b,AOI)); % Angle of incidence loss function
Ee = F1.*((Eb.*F2+ModuleParameters.fd.*Ediff)/E0)*SF; %Effective irradiance
Ee(isnan(Ee))=0; % Set any NaNs to zero

mSAPMResults = pvl_sapm(ModuleParameters, Ee, celltemp);
mSAPMResults.Vmp(mSAPMResults.Vmp<0)=0; %set any spurious negatives to zero
mSAPMResults.Imp(mSAPMResults.Imp<0)=0; %set any spurious negatives to zero
Ms = 28; %Number of modules in series
Mp = 17; %Number of paralell strings
Vdc = Ms *mSAPMResults.Vmp; %Voltage summed by panels in series
Idc = Mp *mSAPMResults.Imp; %Current summed by panels in parallel
Pdc = Vdc .* Idc; %Total array power
Pdc2 = Pdc * 2;
sum_dc = sum(Pdc2); % 1.5546e+07 W of total DC power generation
max_celltemp = max(celltemp);

fig = figure('units','inch','position',[5,5,5,5]);
hold on
plot(doy, celltemp, 'k-')
plot(doy(celltemp == max(celltemp)), max(celltemp), 'ro')
xlim([0,366])
xlabel('Day of the year')
ylabel('Cell Temperature (Degree C)')
title('Cell temperature over the year')
legend('Cell Temperature (Degree C)', append('Maximum Temperature', ...
    num2str(max_celltemp), ' Degree C'), 'Location','NorthWest')
hold off
print(fig,'Q2.png','-dpng','-r800');

%% Question 3
load('SandiaInverterDatabaseSAM2014.1.14.mat')
Inverter = SNLInverterDB(1393);
ACPower = pvl_snlinverter(Inverter,Vdc,Pdc);
ACPower(ACPower<0)=0; %set any spurious negatives to zero
ACPower2 = ACPower * 2;

daily_AC = zeros(366,1);
for i = 1:366
    daily_AC(i) = sum(ACPower2(doy == i));
end

daily_DC = zeros(366,1);
for i = 1:366
    daily_DC(i) = sum(Pdc2(doy == i));
end

DC_max_day = find(daily_DC == max(daily_DC));
DC_min_day = find(daily_DC == min(daily_DC));

fig = figure('units','inch','position',[5,5,5,5]);
hold on
plot(Pdc2(doy == DC_max_day), 'ko-')
plot(ACPower2(doy == DC_max_day), 'bo-')
xlabel('Hour of the day')
ylabel('Power generated (W)')
legend('DC Power Generated', 'AC Power Generated')
title('DC/AC Power generated at max DC power generation day')
hold off
print(fig,'Q3.1.png','-dpng','-r800');

fig = figure('units','inch','position',[5,5,5,5]);
hold on
plot(Pdc2(doy == DC_min_day), 'ko-')
plot(ACPower2(doy == DC_min_day), 'bo-')
xlabel('Hour of the day')
ylabel('Power generated (W)')
legend('DC Power Generated', 'AC Power Generated')
title('DC/AC Power generated at min DC power generation day')
hold off
print(fig,'Q3.2.png','-dpng','-r800');

inverter_effciency = ACPower2 ./ Pdc2;
inverter_effciency(isnan(inverter_effciency))=0;
power_factor = (Inverter.Vac^2 + Pdc.^2).^0.5;
fig = figure('units','inch','position',[5,5,5,5]);
plot(inverter_effciency, power_factor, 'ko', 'MarkerSize',1)
xlabel('Inverter Efficiency')
xlim([0.01, 0.8])
ylabel('Power Factor')
title('Power Factor Vs. Inverter Efficiency')
print(fig,'Q3.3.png','-dpng','-r800');

%% Question 4
it = 10e3:500:500e3; % 10 and 500 kWAC.
ac_sum_list = zeros(length(it),1);
n = 1;
Inverter_temp = Inverter;
for i = it
    multiplier = i/Inverter.Pac0;
    Inverter_temp.Pac0 = Inverter.Pac0 * multiplier;
    Inverter_temp.Pdc0 = Inverter.Pdc0 * multiplier;
    Inverter_temp.Ps0 = Inverter.Ps0 * multiplier;
    Inverter_temp.Pnt = Inverter.Pnt * multiplier;
    ACPower_temp = pvl_snlinverter(Inverter_temp,Vdc,Pdc);
    ACPower_temp(ACPower_temp<0)=0; %set any spurious negatives to zero
    ac_sum_list(n) = sum(ACPower_temp)/10e6;
    n = n + 1;
end

fig = figure('units','inch','position',[5,5,5,5]);
plot(it, ac_sum_list, 'k-')
xlabel('Inverter Capacity (W)')
ylabel('Annual AC power generation (MW)')
print(fig,'Q4.1.png','-dpng','-r800');

n = 1;
pinv = zeros(length(it),1);
for i = it
    pinv(n) = pac2pinv(i);
    n = n + 1;
end

fig = figure('units','inch','position',[5,5,5,5]);
plot(ac_sum_list, pinv, 'b-')
xlabel('Annual AC power generation (MW)')
ylabel('Inverter Price ($/kWAC)')
print(fig,'Q4.2.png','-dpng','-r800');

% There is a huge price drop at Annual AC generation = 0.45 MW,
% Which corresponds to the inverter Pac of approximately 100,000 W
% Thus we can conclude that Xantrex GT 100 inverter is a good fit.

%% Question 5

it = -4:0.1:4;
n = 1;
ac_sum_list = zeros(length(it),1);
for deltaT = it
    celltemp_new = pvl_sapmcelltemp(E, E0, ModuleParameters.a_wind, ...
        ModuleParameters.b_wind, pressure_struct.WindSpeed, ...
        pressure_struct.DryBulb + deltaT, ModuleParameters.delT);
    mSAPMResults = pvl_sapm(ModuleParameters, Ee, celltemp_new);
    mSAPMResults.Vmp(mSAPMResults.Vmp<0)=0; %set any spurious negatives to zero
    mSAPMResults.Imp(mSAPMResults.Imp<0)=0; %set any spurious negatives to zero
    Ms = 28; %Number of modules in series
    Mp = 17; %Number of paralell strings
    Vdc = Ms *mSAPMResults.Vmp; %Voltage summed by panels in series
    Idc = Mp *mSAPMResults.Imp; %Current summed by panels in parallel
    Pdc = Vdc .* Idc; %Total array power
    Inverter = SNLInverterDB(1393);
    ACPower = pvl_snlinverter(Inverter,Vdc,Pdc);
    ACPower(ACPower<0)=0; %set any spurious negatives to zero
    ac_sum_list(n) = sum(ACPower);
    n = n + 1;
end

fig = figure('units','inch','position',[5,5,5,5]);
hold on
plot(it, ac_sum_list./1e6, 'k-')
plot(-4, ac_sum_list(1)/1e6, 'ro')
plot(-2, ac_sum_list(21)/1e6, 'ro')
plot(-0, ac_sum_list(41)/1e6, 'ro')
plot(2, ac_sum_list(61)/1e6, 'ro')
plot(4, ac_sum_list(81)/1e6, 'ro')
xlabel('Ambient Temperature Change (Degree C)')
ylabel('Annual Power Generation (MW)')
legend('Trend', 'Temperature Change')
hold off
print(fig,'Q5.png','-dpng','-r800');

%% Question 6
real_struct = load('realACPower2016.mat');
ACPower_real = real_struct.PowerAC;
ACPower_real(ACPower_real<0) = 0;

daily_AC_real = zeros(366,1);
for i = 1:366
    daily_AC_real(i) = sum(ACPower_real(doy == i));
end

fig = figure('units','inch','position',[5,5,12,5]);
subplot(1,2,1)
hold on
plot(daily_AC_real, 'k-')
plot(daily_AC, 'b-')
legend('Measured AC Power Generation', 'Modeled AC Power Generation')
title('Original')
hold off

subplot(1,2,2)
hold on
plot(daily_AC_real * 30, 'k-')
plot(daily_AC, 'b-')
legend('Measured AC Power Generation', 'Modeled AC Power Generation')
title('30x')
hold off
print(fig,'Q6.1.png','-dpng','-r800');

%% Done
close all;



