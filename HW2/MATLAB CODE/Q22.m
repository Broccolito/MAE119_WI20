clc; clear;
Y = 12 * cosd(10);
x = 180 - atand(Y/15);
Phi = 33;

Wst_mat = [];
for d = 1:30
    Wst_list = [];
    for i = 1:12
        today_time = datetime(2019,i,d,0:23,0,0);
        Time = pvl_maketimestruct(datenum(today_time), -8);
        Location = pvl_makelocationstruct(Phi, -117);
        [SunAz, SunFl, AppSunFl, SolarTime] = pvl_ephemeris(Time, Location);
        WSt = interp1(SunAz, SolarTime, x);
        Wst_list = [Wst_list, WSt];
    end
    Wst_mat = [Wst_mat; Wst_list];
end

tlist = [];
for i = 1:12
    
   tlist = [tlist; Wst_mat(:,i)];
    
end


plot(1:360, tlist)
title('Waking up time over the year at San Diego')
xlabel('Day')
ylabel('Wake Up Time')

