%% p.2
clear all;
clc;

% at z=20m
u=8;
sd=3.7;

syms c k

ksim=vpasolve(sd^2==(u^2)*((gamma(1+2/k)/gamma(1+1/k)^2)-1),k);
ksim=double(ksim);
numc=vpasolve(u==c*gamma(1+1/ksim),c);
numc=double(numc);

k=ksim;
c=numc;
%Check
[m,v]=wblstat(c,k);

disp('At z = 20 m, the values of k and c are: ')
disp(k)
disp(c)

% at z=100m
u=6.357;
sd=2.94;

syms c k

ksim=vpasolve(sd^2==(u^2)*((gamma(1+2/k)/gamma(1+1/k)^2)-1),k);
ksim=double(ksim);
numc=vpasolve(u==c*gamma(1+1/ksim),c);
numc=double(numc);

k=ksim;
c=numc;
%Check
[m,v]=wblstat(c,k);

disp('At z = 100 m, the values of k and c are: ')
disp(k)
disp(c)

%% p.3
clc;
clear;
u2=[3.2 4.5 2.9 5.2 6.1 1.5 4.2]; %z=10m
u1=[4 5.5 3.6 6.4 7.5 1.8 5.2]; %z=30m

%check alpha
a=.16669;
x=1:.5:max(u1);
y = x.*3^a;

%compare results from table and alpha
fig = figure('units','inch','position',[5,5,6,5]);
hold on
scatter(u2,u1, 'MarkerFaceColor', 'k'); 
plot(x,y);
title('Wind Speeds at 10m and 30m')
xlabel('Wind speed at 10m (m/s)')
ylabel('Wind speed at 30m (m/s)')
legend('Windspeed Measured', 'Windspeed Predicted')
hold off
print(fig,'fig3.png','-dpng','-r800');