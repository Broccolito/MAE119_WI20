%% Import pvlib
addpath(genpath('PVLib 1.4 Release'));
%% Setup location and time to analyze
today_time=datetime(2019,6,20,0:23,0,0);
% Feed in time and time zone
Time = pvl_maketimestruct(datenum(today_time),-8);
Location = pvl_makelocationstruct(32.86,-117.22); %San Diego lat and lon
%% Obtain sun position angles
[SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location);
SunZen = 90 - SunEl;

%% Obtain air mass
AMa = pvl_relativeairmass(SunZen);
%there is also an absolute air mass function, requires pressure

%% Plot air mass
plot(today_time,AMa); xlabel('Time'); ylabel('Air mass');

