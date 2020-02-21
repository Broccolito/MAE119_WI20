%% Import pvlib
addpath(genpath('PVLib 1.4 Release'));
%% Setup location and time to analyze
today_time=datetime(2019,12,21,0:23,0,0);
% Feed in time and time zone
Time = pvl_maketimestruct(datenum(today_time),-8);
Location = pvl_makelocationstruct(32.86,-117.22); %San Diego lat and lon
%% Obtain sun position angles
[SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location);
SunZen=90 - SunEl;
%% Plot Solar angles
plot(today_time,SunAz,today_time,SunEl,today_time,SunZen); grid on
legend('Azimuth angle','Elevation angle','Zenith angle')
xlabel('Hour of day')
ylabel('Angle (deg)')
