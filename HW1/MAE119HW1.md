Problem 1

4. Using the MATLAB code below:

   ```matlab
   % Q1.4
   
   disp('San Diego')
   n = 359;
   delta = 23.45 * sind(360*(284+n)/365);
   omega = -1.5*15;
   lat = 33;
   zenSD = acosd(cosd(lat)*cosd(delta)*cosd(omega)+sind(lat)*sind(delta));
   betaa = 0;
   betab = 30;
   betac = 45;
   azib = 180;
   azic = 150;
   
   incia = acosd(sind(delta)*sind(lat)*cosd(betaa)-...
       sind(delta)*cosd(lat)*sind(betaa)*cosd(azib)+...
       cosd(delta)*cosd(lat)*cosd(betaa)*cosd(omega)+...
       cosd(delta)*sind(lat)*sind(betaa)*cosd(azib)*cosd(omega)+...
       cosd(delta)*sind(betaa)*sind(azib)*sind(omega));
   
   incib = acosd(sind(delta)*sind(lat)*cosd(betab)-...
       sind(delta)*cosd(lat)*sind(betab)*cosd(azib)+...
       cosd(delta)*cosd(lat)*cosd(betab)*cosd(omega)+...
       cosd(delta)*sind(lat)*sind(betab)*cosd(azib)*cosd(omega)+...
       cosd(delta)*sind(betab)*sind(azib)*sind(omega));
   
   incic = acosd(sind(delta)*sind(lat)*cosd(betac)-...
       sind(delta)*cosd(lat)*sind(betac)*cosd(azib)+...
       cosd(delta)*cosd(lat)*cosd(betac)*cosd(omega)+...
       cosd(delta)*sind(lat)*sind(betac)*cosd(azib)*cosd(omega)+...
       cosd(delta)*sind(betac)*sind(azib)*sind(omega));
   
   disp(incia)
   disp(incib)
   disp(incic)
   
   
   disp('Melbourne')
   n = 359;
   delta = 23.45 * sind(360*(284+n)/365);
   omega = -1.5*15;
   lat = 37;
   zenSD = acosd(cosd(lat)*cosd(delta)*cosd(omega)+sind(lat)*sind(delta));
   betaa = 0;
   betab = 30;
   betac = 45;
   azib = 180;
   azic = 150;
   
   incia = acosd(sind(delta)*sind(lat)*cosd(betaa)-...
       sind(delta)*cosd(lat)*sind(betaa)*cosd(azib)+...
       cosd(delta)*cosd(lat)*cosd(betaa)*cosd(omega)+...
       cosd(delta)*sind(lat)*sind(betaa)*cosd(azib)*cosd(omega)+...
       cosd(delta)*sind(betaa)*sind(azib)*sind(omega));
   
   incib = acosd(sind(delta)*sind(lat)*cosd(betab)-...
       sind(delta)*cosd(lat)*sind(betab)*cosd(azib)+...
       cosd(delta)*cosd(lat)*cosd(betab)*cosd(omega)+...
       cosd(delta)*sind(lat)*sind(betab)*cosd(azib)*cosd(omega)+...
       cosd(delta)*sind(betab)*sind(azib)*sind(omega));
   
   incic = acosd(sind(delta)*sind(lat)*cosd(betac)-...
       sind(delta)*cosd(lat)*sind(betac)*cosd(azib)+...
       cosd(delta)*cosd(lat)*cosd(betac)*cosd(omega)+...
       cosd(delta)*sind(lat)*sind(betac)*cosd(azib)*cosd(omega)+...
       cosd(delta)*sind(betac)*sind(azib)*sind(omega));
   
   disp(incia)
   disp(incib)
   disp(incic)
   ```

   The irradiance angles are found to be    60.3315, 88.2067 and 102.2376 degrees at San Diego and

      64.0022, 91.9518, 105.9659 degrees in Melbourne.



Problem 2:

1. For the place where I am currently at, the latitude is 32.860 N, and the longitude is 117.220 W.

2. By executing the code:

   ```matlab
   %% Import pvlib
   addpath(genpath('PVLib 1.4 Release'));
   %% Setup location and time to analyze
   today_time=datetime(2019,6,20,0:23,0,0);
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
   ```

   We are able to plot the angles at June 20th (summer solstice)

   ![image-20200122114058450](C:\Users\wanju\AppData\Roaming\Typora\typora-user-images\image-20200122114058450.png)

   By executing the code:

   ```matlab
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
   ```

   We are able to plot the angles at Dec 21st (winter solstice)

   ![image-20200122114205399](C:\Users\wanju\AppData\Roaming\Typora\typora-user-images\image-20200122114205399.png)



3. By executing the code:

   ```matlab
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
   ```

   and

   ```matlab
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
   
   %% Obtain air mass
   AMa = pvl_relativeairmass(SunZen);
   %there is also an absolute air mass function, requires pressure
   
   %% Plot air mass
   plot(today_time,AMa); xlabel('Time'); ylabel('Air mass');
   ```

   We are able to plot the air mass as a function of the time of the day of the summer solstice and winter solstice. 

   As we can see from the graph, the air mass is smaller during summer times. This is due to the fact that sunlight goes through the atmosphere more directly during the summer whereas less directly during the winter.

   ![image-20200122114600480](C:\Users\wanju\AppData\Roaming\Typora\typora-user-images\image-20200122114600480.png)

   ![image-20200122114631191](C:\Users\wanju\AppData\Roaming\Typora\typora-user-images\image-20200122114631191.png)

   

   

   4. By executing the code:

      ```matlab
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
      ```

      We are able to plot the relation between angle of incidence and time of the day with different solar panel configurations.

      ![image-20200122115015776](C:\Users\wanju\AppData\Roaming\Typora\typora-user-images\image-20200122115015776.png)

      ![image-20200122115026972](C:\Users\wanju\AppData\Roaming\Typora\typora-user-images\image-20200122115026972.png)

      ![image-20200122115040230](C:\Users\wanju\AppData\Roaming\Typora\typora-user-images\image-20200122115040230.png)

      As we can clearly see, the first configuration (20 degree tilt, south) is better than the second configuration (30 degree tilt, west). However, the configuration that I propose (37 degree tilt, south) should have better outcomes in the winter than both of the previous configurations due to a smaller angle of incidence.

   

