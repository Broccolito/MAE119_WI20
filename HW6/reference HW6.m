%% HW 6
%Reuben Lumaban
%Eni Ikuku


%% p.2
clear all;
clc;
% at z=20m
% u=8;
% sd=3.7;

% at z=100m
u=6.357;
sd=2.94;

syms c k

% u=c*gamma(1+1/k);
% sd^2=(u^2)*((gamma(1+2/k)/gamma(1+1/k)^2)-1);

numk=vpasolve(sd^2==(u^2)*((gamma(1+2/k)/gamma(1+1/k)^2)-1),k);
numk=double(numk);
numc=vpasolve(u==c*gamma(1+1/numk),c);
numc=double(numc);

k=numk;
c=numc;
%Check
[m,v]=wblstat(c,k);

%% p.3
clc;
clear;
u2=[3.2 4.5 2.9 5.2 6.1 1.5 4.2]; %z=10m
u1=[4 5.5 3.6 6.4 7.5 1.8 5.2]; %z=30m

%check alpha
a=.16669;
%y=1:.5:max(u1);
x=1:.5:max(u1);
for i=1:length(x)
    y(i)=x(i)*3^a;
end

%compare results from table and alpha
figure 
hold on
scatter(u2,u1); 
plot(x,y);
title('Wind Speeds at different heights')
xlabel('Wind speed at 10m (m/s)')
ylabel('Wind speed at 30m (m/s)')
legend('Measured windspeed', 'Predicted windspeed using alpha parameter')