clear;
close all;
warning off
addpath(genpath('PVLib'))

power_load = readtable('Konocti_Load.csv');
power_load = table2array(power_load(:,3));

calc_power_generation

solar_range = 10:80; %solar capacity, MW
battery_range = 10:200; %battery capacity, MWh
life_cycle_cost = zeros(200,80);

for solar = solar_range
for battery = battery_range

total_power_generation = power_generation*solar; %solar generation * number of arrays, MW

Battery = zeros(288,1); %record charge of battery each hour
Battery(1) = battery;
PowerOut = zeros(288,1); %record number of unserved MW (will sum later)

for hour = 1:length(power_load)
if total_power_generation(hour) < power_load(hour)
Battery(hour+1) = Battery(hour) - (power_load(hour)-total_power_generation(hour)); %charge battery
if Battery(hour+1) < 0
    PowerOut(hour) = abs(Battery(hour+1)); %record # of unserved MWh
    Battery(hour+1) = 0; %reset battery to zero
end

elseif total_power_generation(hour) >= power_load(hour)
Battery(hour+1) = Battery(hour) + (total_power_generation(hour) - power_load(hour));
if Battery(hour+1) > battery
    Battery(hour+1) = battery;
end
end
end



power_outage_cost = zeros(21,1); 
investment = zeros(21,1); 
OnM = zeros(21,1); 
incentives = zeros(21,1); 

power_outage_cost(1) = sum(PowerOut)*1000*14;
investment(1) = 1800*solar*1000 + 600*battery*1000; 
investment(11) = 600*battery*1000; 
OnM(2) = 9*solar*1000; 
incentives(1) = -1800*solar*1000*0.3;

for m = 3:21
    OnM(m) = OnM(2)*(1.05)^(m-2);
end

for p = 2:21
    power_outage_cost(p) = power_outage_cost(1)*(1.05)^(p-2);
end

sumcost = investment + power_outage_cost + OnM + incentives;

discount = zeros(21,1);
for n = 1:21
    discount(n) = sumcost(n)/(1.06)^(n-1);
end

life_cycle_cost(battery,solar) = sum(discount);
end
end

life_cycle_cost = life_cycle_cost(min(solar_range):end, ...
    min(battery_range):end);

figure('units','inch','position',[5,5,6,5]);
surf(solar_range, battery_range, life_cycle_cost)
% print(fig,'fig2.png','-dpng','-r800');

fig = figure('units','inch','position',[5,5,6,5]);
hold on
contourf(solar_range,battery_range,...
    life_cycle_cost,'LineStyle','none')
colormap(summer)
xlabel('Solar Panel Capacity')
ylabel('Battery Capacity')
hold off
print(fig,'fig2.png','-dpng','-r800');

minimum = min(min(life_cycle_cost));
[b,s] = find(life_cycle_cost == minimum);

battery_best = battery_range(b);
solar_best = solar_range(s);

disp('Best Battery Capacity:  ')
disp(battery_best);
disp('Best Solar Panel Capacity:  ')
disp(solar_best);
