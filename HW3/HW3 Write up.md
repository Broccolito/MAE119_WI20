##### Wanjun Gu

##### A13487962

##### MAE119

##### HW3



#### Question 1

The plots can be generated using the code:

```matlab
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
```

![Q1.1](C:\Users\wanju\Desktop\MAE 119\HW3\Q1.1.png)

![Q1.2](C:\Users\wanju\Desktop\MAE 119\HW3\Q1.2.png)

![Q1.3](C:\Users\wanju\Desktop\MAE 119\HW3\Q1.3.png)



#### Question 2

Using the code:

```matlab
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
```

We get that the total DC power generated is 1.5546e+07 W.

The cell temperature over the year is graphed as:

![Q2](C:\Users\wanju\Desktop\MAE 119\HW3\Q2.png)



#### Question 3

Using the code:

```matlab
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
```



We know that the days with the highest power generation is the 143rd day of the year. The day with the smallest power generation is the fifth day of the year.

At the day with maximum power generation:

![Q3.1](C:\Users\wanju\Desktop\MAE 119\HW3\Q3.1.png)

At the day with minimal power generation:

![Q3.2](C:\Users\wanju\Desktop\MAE 119\HW3\Q3.2.png)

The power factor Vs. convertor efficiency and be shown as:

![Q3.3](C:\Users\wanju\Desktop\MAE 119\HW3\Q3.3.png)



#### Question 4

Using the code:

```matlab
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
```

We know that the relation between total annual AC power generation in MWh/a versus the temperature change can be plotted as:

![Q4.1](C:\Users\wanju\Desktop\MAE 119\HW3\Q4.1.png)

And the cost can be plotted as:

![Q4.2](C:\Users\wanju\Desktop\MAE 119\HW3\Q4.2.png)

As the plots show, There is a huge price drop at Annual AC generation = 0.45 MW, which corresponds to the inverter Pac of approximately 100,000 W. Thus we can conclude that Xantrex GT 100 inverter is a good fit.



#### Question 5

Using the code:

```matlab
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

slope = range(ac_sum_list) / range(it);
% The panel temperature coefficeint remains constant
% over different temperature ranges. The coefficient
% is indicated by the slope of plot 5.1, which is
% 2.7489e+04 W/degree C. This means that by one deg
% temperature increase, the cell will tend to generate
% 2.7489e+04 W less power a year.
```

We can plot:

![Q5](C:\Users\wanju\Desktop\MAE 119\HW3\Q5.png)

On this graph we can see that there is clearly a negative correlation between annual power generation and the ambient temperature change. The panel temperature coefficient remains constant over different temperature ranges. The coefficient is indicated by the slope of plot 5.1, which is 2.7489e+04 W/degree C. This means that by one degree temperature increase, the cell will tend to generate 2.7489e+04 W less power a year.



#### Question 6

Using the code:

```matlab
%% Question 6
real_struct = load('realACPower2016.mat');
ACPower_real = real_struct.PowerAC;
ACPower_real(ACPower_real<0) = 0;

daily_AC_real = zeros(366,1);
for i = 1:366
    daily_AC_real(i) = sum(ACPower_real(doy == i));
end
daily_AC_real = daily_AC_real * 28;

fig = figure('units','inch','position',[5,5,12,5]);
hold on
plot(daily_AC_real./28, 'c-')
plot(daily_AC_real, 'k-')
plot(daily_AC, 'b-')
legend('Single-array Measured AC Power Generation',...
    'Total Measured AC Power Generation',...
    'Modeled AC Power Generation')
xlabel('Day of the year')
ylabel('AC Power Generation (W)')
hold off
print(fig,'Q6.1.png','-dpng','-r800');

array_yield = sum(daily_AC_real./28);
final_yield = sum(daily_AC_real);
reference_yield = sum(daily_AC);
panel_efficiency_measure = final_yield / (sum(AOI)*1000);
panel_efficiency_model = reference_yield / (sum(AOI)*1000);
% inverter_efficiency is not obtainable in real
% measurements because there is no data on DC 
% power output
performance_ratio = final_yield / reference_yield;
capacity_factor_measure = final_yield ./ (max(daily_AC_real) .* 366);
capacity_factor_model = reference_yield ./ (max(daily_AC) .* 366);
```

We can first show the graph of the real measured power generation Vs. the modeled one:

![Q6.1](C:\Users\wanju\Desktop\MAE 119\HW3\Q6.1.png)

We can also easily calculate that:

The array yield of the measured panel is 2.9983e+05 W

The final yield of the measured panel is 8.3953e+06 W

The reference yield of the model panel is 8.9046e+06 W

The panel efficiency of the measured one is 1.06%

The panel efficiency of the model is 1.13%

The capacity factor for the measured panel is 0.6185

The capacity factor for the model is 0.5865

There is however no way to calculate the inverter profile with the measured data given because the DC output is not given. Also, assumptions are made when processing the measured data that the AOI used in the model is the same as real life.



As we can see, the model data and measure data match well. However, in most cases, the model slightly overestimate the power generation. This is likely caused by weather or dust on the panel. The efficiency of the panel is low. This may due to the fact that the tilting angle is not optimized for the latitude. However, more reasons behind low efficiency need to be further elucidated.