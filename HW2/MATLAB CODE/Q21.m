clc; clear;
Y = 12 * cosd(10);
x = 180 - atand(Y/15);
Phi = 33;
today_time = datetime(2019,2,1,0:23,0,0);
Time = pvl_maketimestruct(datenum(today_time), -8);
Location = pvl_makelocationstruct(Phi, -117);
[SunAz, SunFl, AppSunFl, SolarTime] = pvl_ephemeris(Time, Location);
WSt = interp1(SunAz, SolarTime, x);


% WSt = 9.7177
% Sun hits window at 9:43 a.m