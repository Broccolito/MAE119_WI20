%% Import pvlib
addpath(genpath('PVLib 1.4 Release'));
%% Setup location and time to analyze
today_time=datetime(2019,12,21,0:23,0,0);
% Feed in time and time zone
Time = pvl_maketimestruct(datenum(today_time),-8);
Location = pvl_makelocationstruct(32.86,-117.22); %San Diego lat and lon
%% Obtain sun position angles
[SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location);
SunZen = 90 - SunEl;

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup solar panel orientation and get angle of incidence
SurfTilt1 = 20; % Array tilt angle (deg)
SurfAz1 = 180; %Array azimuth (180 deg indicates array faces South)
AOI1 = pvl_getaoi(SurfTilt1,SurfAz1,SunZen,SunAz);
%% Plot AOI
plot(today_time, AOI1); xlabel('Time'); ylabel('Angle of incidence (deg)')
saveas(gcf, '20_south.png')

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup solar panel orientation and get angle of incidence
SurfTilt2 = 30; % Array tilt angle (deg)
SurfAz2 = 270; %Array azimuth (180 deg indicates array faces South)
AOI2 = pvl_getaoi(SurfTilt2,SurfAz2,SunZen,SunAz);
%% Plot AOI
plot(today_time, AOI2); xlabel('Time'); ylabel('Angle of incidence (deg)')
saveas(gcf, '30_west.png')

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup solar panel orientation and get angle of incidence
SurfTilt3 = 37; % Array tilt angle (deg)
SurfAz3 = 180; %true south
AOI3 = pvl_getaoi(SurfTilt3,SurfAz3,SunZen,SunAz);
%% Plot AOI
plot(today_time, AOI3); xlabel('Time'); ylabel('Angle of incidence (deg)')
saveas(gcf, '37_south.png')
