addpath(genpath('PVLib 1.4 Release'))

month = 1:1:12;
for imonth = 1:12
    today_time = datetime(2019, imonth, 01, 0:23, 0, 0);
    % Time and time zone
    Time = pvl_maketimestruct(datenum(today_time),-8);
    Location = pvl_makelocationstruct(33, -117);
    
    [SunAz, SunEl,AppSunEl,SolarTime] = pvl_ephemeris(Time, Location);
    sun_azimuth(:,imonth) = SunAz;
    solar_time(:,imonth) = SolarTime;
end

% Plot solar position
plot(sun_azimuth, solar_time)
legend(month)
hold on
xlim([120 180])