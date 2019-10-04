clc
clear
close all

%% Main procedure to generate measurements

time=200;%simulation time
deltat=0.05;%measurement sample time

%initial condition for the system
xhat0=[6500.4,349.14,-1.8093,-6.7967,.6932]'; %mean of the initial condition
phat0=diag([1e-6,1e-6,1e-6,1e-6,0]); %covariance matrix of the initial condition
x0= mvnrnd(xhat0,phat0,1); %initial condition vector
xsim=x0;

%solving system state equations
options=odeset('AbsTol',1e-7,'RelTol',1e-7);
Q=diag([0,0,2.4064e-5,2.4064e-5,1e-6])*0.05; %covariance matrix of noise
for i=1:(time/deltat)
    [T,X]=ode45(@fdot4,[i*deltat deltat*(i+1)],xsim,options);
    w= mvnrnd([0 0 0 0 0]',Q,1);
    xsim= mvnrnd(X(end,1:5),Q,1);
    range=((xsim(1)-6374)^2+xsim(2)^2)^0.5;
    bangle=atan((xsim(2)-0)/(xsim(1)-6374));
    v = mvnrnd([0 0]',[0.001 0;0 0.017],1);
    state(:,i)=(X(end,1:5))';
    y(i,1)=range+v(1);
    y(i,2)=bangle+v(2);
end

