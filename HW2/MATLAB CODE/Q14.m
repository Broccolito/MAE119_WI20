%% Add pvlib to path
addpath(genpath('pvlib'));
%% Setup location and time to analyze
today_time=datetime(2019,1:12,7,11,0,0);
Time = pvl_maketimestruct(datenum(today_time),-8);
Location = pvl_makelocationstruct(33,-117,0);
%% How can we get when is solar noon?
[SunAz, SunEl, AppSunEl, SolarTime] = pvl_ephemeris(Time,Location);
[~,isolnoon]=min(abs(SolarTime-12)); % find when solartime is closer to 12
tsolnoon=today_time(isolnoon);
Time = pvl_maketimestruct(datenum(tsolnoon),-8); %redefine time to be solar noon
%% Get clear sky GHI using Ineichen-Perez clear sky model
[ClearSkyGHI_1, ClearSkyDNI_1, ClearSkyDIF_1]= pvl_clearsky_ineichen(Time,Location,1);
% [ClearSkyGHI_7, ClearSkyDNI_7, ClearSkyDIF_7]= pvl_clearsky_ineichen(Time,Location,7);
% fprintf(['GHI for TL=1 is ',num2str(ClearSkyGHI_1),', kt=',num2str(840/ClearSkyGHI_1),'\n'])
% fprintf(['GHI for TL=7 is ',num2str(ClearSkyGHI_7),', kt=',num2str(840/ClearSkyGHI_7),'\n'])
% %How should DNI and DHI change with TL?
% fprintf(['DNI for TL=1 is ',num2str(ClearSkyDNI_1),'\n'])
% fprintf(['DNI for TL=7 is ',num2str(ClearSkyDNI_7),'\n'])
% fprintf(['DIF for TL=1 is ',num2str(ClearSkyDIF_1),'\n'])
% fprintf(['DIF for TL=7 is ',num2str(ClearSkyDIF_7),'\n'])
ClearSkyGHI_1

%% GHI tends to increase with respect to increasing altitude and it tends to decrease with increasing T_l or airmass