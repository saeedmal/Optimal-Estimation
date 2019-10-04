
%% Main procedure to generate measurements
clc
clear
close all

global R0 rv muk

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

save('measurement.mat','state','y');

disp('Measurement generation done');

%% Extended Kalman Filter(EKF)
load measurement %load measurement information
 
jk=2;
deltat=0.05*jk; 
rind=jk; %index of multiplication
for iii=1:(T/deltat) %sampling from the measurement to go through filters
    state1(:,iii)=state(:,rind*iii);
    y1(iii,:)=y(rind*iii,:);
end

time=200; %simulation time

%initial condition for the system
xhat0=[6500.4 349.14 -1.8093 -6.7967 0]; %mean of the initial condition
phat0=diag([1e-6,1e-6,1e-6,1e-6,1]); %covariance matrix of the initial condition
x0=[xhat0 (phat0(:,1))' (phat0(:,2))' (phat0(:,3))' (phat0(:,4))' (phat0(:,5))'];%initial state1

%solving system state equations
options=odeset('AbsTol',1e-7,'RelTol',1e-7);
muk=xhat0;
pk=phat0;
xsim=x0;
for i=1:(time/deltat)

    [T,X]=ode45(@fdot5,[i*deltat deltat*(i+1)],xsim,options);

    muk=(X(end,1:5))';
    pk=[(X(end,6:10))' (X(end,11:15))' (X(end,16:20))' (X(end,21:25))' (X(end,26:30))'];

    %y=h(x)+v,H=dh/dx linearization of the state space
    h1=((muk(1) - R0)^2 + muk(2)^2)^(1/2);
    h2=atan(muk(2)/(muk(1) - R0));
    H(1,1)=(muk(1) - R0)/h1; 
    H(1,2)=muk(2)/h1;
    H(1,3:5)=0;
    H(2,1)=-muk(2)/(h1^2);
    H(2,2)=(muk(1) - R0)/(h1^2);
    H(2,3:5)=0;

    %calculating kalman gain
    K=pk*H'*inv((H*pk*H')+rv);
    muk=muk+K*((y1(i,1:2))'-[h1;h2]);
    pk=pk-K*H*pk;

    %update for xsim
    xsim=[muk' (pk(:,1))' (pk(:,2))' (pk(:,3))' (pk(:,4))' (pk(:,5))'];

    xe(:,i)=muk; %saving mean for the current j
    Pe(:,i)=diag(pk);   %saving covariance for the current j

end


figure
t=[deltat:deltat:time]; %time vector
real_measurement=state(:,rind:rind:end);
plot(t,xe(1,:)-real_measurement(1,:),'linewidth',2);
hold on
plot(t,3*sqrt(Pe(rind:rind:end,:)),'r');%,t,-3*sqrt(Pe(rind:rind:end,:)),'r','linewidth',3);
rmse=rms(xe(1,:)-real_measurement(1,:));
title(['error for the 1st state, RMSE= ',num2str(rmse)]);
grid on;

disp('Done with EKF')
keyboard
%% This section implements UKF 
for jk=1:8
    
    load measurement
    
    dt=0.05;
    dtm=dt*jk;
    Q=diag([2.4064e-5,2.4064e-5,1e-6])*dtm;%process noise Q
    
    L=5;        %number of states                        
    m=2;        %size of each measurement vector  
    
    %scaling factors                                 
    ki=0;                                                                           
    s=zeros(L,2*L+1);
    rv=[0.001 0;0 0.017];
    musim=xhat0;
    psim=phat0;
    
    options=odeset('AbsTol',1e-7,'RelTol',1e-7);
    for i=1:(time/dtm)
    
        [u,uw]=sigmap(musim,psim);%xsim and psim will be updated one time for each i
        
        for  j=1:(2*L+1) 
            [T,Xk]=ode45(@fdot4,[i*dtm dtm*(i+1)],u(:,j),options);
            s(:,j)=(Xk(end,:))';%save sigma points for each i:x(k+1|k)
        end
        
        %generating Wi for mean and covariance
        muk=zeros(L,1);yhatk=zeros(m,1);
        pxx=zeros(L,L);
        pxy=zeros(L,m);
        pyy=zeros(m,m);
    
        for k=1:(2*L+1)
                muk=muk+uw(k)*s(:,k);
                yhatk=yhatk+uw(k)*out(s(:,k));
        end
        
        for kk=1:(2*L+1)
          pxx=pxx+uw(kk)*(s(:,kk)-muk)*(s(:,kk)-muk)';
          pxy=pxy+uw(kk)*(s(:,kk)-muk)*(out(s(:,kk))-yhatk)';
          pyy=pyy+uw(kk)*(out(s(:,kk))-yhatk)*(out(s(:,kk))-yhatk)';
        end
        
        pyy=pyy+rv;
        
        gg=[zeros(1,3);zeros(1,3);eye(3)];
        gQg=gg*Q*gg';
        pxx=pxx+gQg; 
        %updating mukplus1 and sigkplus1 using measurment
        K=pxy*inv(pyy);%kalman gain
        rind=round(dtm/dt);
        yk=(y(rind*i,:))';%measurement
        musim=muk+K*(yk-yhatk); %update for mean using measurement musim=mu(k+1|k+1)
        psim=pxx-K*pxy'; %updating covariance using measurement
        xe(:,i)=musim; %save mean
        Pe(:,i)=diag(psim); %save covariance
    end
end

%% This section is the algorithm for CUT4
musim=xhat0;
psim=phat0;
options=odeset('AbsTol',1e-7,'RelTol',1e-7);
for i=1:(time/dtm)
    [u,uw]=sigmap(musim,psim);%xsim and psim will be updated one time for each i
    
    n=size(u,2);
    for  j=1:n
        [T,Xk]=ode45(@fdot4,[i*dtm dtm*(i+1)],u(:,j),options);
        s(:,j)=(Xk(end,:))';%save sigma points for each i:x(k+1|k)
    end
    %generating Wi for mean and covariance
    muk=zeros(L,1);yhatk=zeros(m,1);
    pxx=zeros(L,L);
    pxy=zeros(L,m);
    pyy=zeros(m,m);
    for k=1:(n)
        muk=muk+uw(k)*s(:,k);
        yhatk=yhatk+uw(k)*out(s(:,k));
    end
    for kk=1:(n)
        pxx=pxx+uw(kk)*(s(:,kk)-muk)*(s(:,kk)-muk)';
        pxy=pxy+uw(kk)*(s(:,kk)-muk)*(out(s(:,kk))-yhatk)';
        pyy=pyy+uw(kk)*(out(s(:,kk))-yhatk)*(out(s(:,kk))-yhatk)';
    end
    
    pyy=pyy+rv;
    gg=[zeros(1,3);zeros(1,3);eye(3)];
    gQg=gg*Q*gg';
    
    pxx=pxx+gQg; %propagate
    
    %updating mukplus1 and sigkplus1 using measurment
    K=pxy*inv(pyy);%kalman gain
    rind=round(dtm/dt);
    yk=(y(rind*i,:))';%measurement
    musim=muk+K*(yk-yhatk); %update for mean using measurement musim=mu(k+1|k+1)
    psim=pxx-K*pxy'; %updating covariance using measurement
    
    xe(:,i)=musim; %save mean
    Pe(:,i)=diag(psim); %save covariance
end

t=[dtm:dtm:time]; %time vector

state=state(:,rind:rind:end);

subplot(3,2,1);plot(t,xe(1,:)-state(1,:),'linewidth',2);
hold on
plot(t,3*sqrt(Pe(1,:)),'r',t,-3*sqrt(Pe(1,:)),'r','linewidth',3);
title('error for the 1st state');
rmse(1)=rms(xe(1,:)-state(1,:));
grid on;

%re-entry model state space
function xdot=fdot4(t,x)
%model parameters
beta0=0.59783;
R0=6374;
Gm0=3.986e5;
H0=13.406;
R=(x(1)^2+x(2)^2)^0.5;
G=-Gm0/(R^3);
V=(x(3)^2+x(4)^2)^0.5;
beta=beta0*exp(x(5));
D=-beta*V*exp((R0-R)/H0);
%model state space
xdot(1,1)=x(3);
xdot(2,1)=x(4);
xdot(3,1)=D*x(3)+G*x(1);
xdot(4,1)=D*x(4)+G*x(2);
xdot(5,1)=0;
end

%re-entry model state space
function xdot=fdot5(t,x)
global R0 rv muk
%model parameters
beta0=0.59783;
R0=6374;
rv=[0.001 0;0 0.017];
Gm0=3.986e5;
H0=13.406;
R=(x(1)^2+x(2)^2)^0.5;
G=-Gm0/(R^3);
V=(x(3)^2+x(4)^2)^0.5;
beta=beta0*exp(x(5));
D=-beta*V*exp((R0-R)/H0);
%model state space
%partI: mudot=f(mu)
xdot(1,1)=x(3);
xdot(2,1)=x(4);
xdot(3,1)=D*x(3)+G*x(1);
xdot(4,1)=D*x(4)+G*x(2);
xdot(5,1)=0;
%part II: sigmadot=A*sigma+sigma'*A+Q
A=zeros(5,5);
%evaluate A[5*5]=df/dx at muk
R1=(muk(1)^2+muk(2)^2)^0.5;
G1=-Gm0/(R1^3);
V1=(muk(3)^2+muk(4)^2)^0.5;
D1=-beta*V1*exp((R0-R1)/H0);
diffD=[-D*muk(1)/(H0*R1);
 (beta0*muk(2)*exp(muk(5))*exp((R0 - R1)/H0)*V1)/(H0*R);
  -(beta0*muk(3)*exp(muk(5))*exp((R0 - R1)/H0))/V1;
  -(beta0*muk(4)*exp(muk(5))*exp((R0 - R1)/H0))/V1;
  D1];
diffG=[(3*Gm0*muk(1))/R1^5;
       (3*Gm0*muk(2))/R1^5;
                      0;
                      0;
                      0];
A(1,:)=[0 0 1 0 0];
A(2,:)=[0 0 0 1 0];
A(5,:)=[0 0 0 0 0];
%partial diff of D
A(3,:)=[diffD(1)*muk(3)+diffG(1)*muk(1)+G1;
        diffD(2)*muk(3)+diffG(2)*muk(1);
        diffD(3)*muk(3)+diffG(3)*muk(1)+D1;
        diffD(4)*muk(3)+diffG(4)*muk(1);
        diffD(5)*muk(3)+diffG(5)*muk(1)]';
A(4,:)=[diffD(1)*muk(4)+diffG(1)*muk(2);
        diffD(2)*muk(4)+diffG(2)*muk(2)+G1;
        diffD(3)*muk(4)+diffG(3)*muk(2);
        diffD(4)*muk(4)+diffG(4)*muk(2)+D1;
        diffD(5)*muk(4)+diffG(5)*muk(2)]';
Q=diag([2.4064e-5,2.4064e-5,1e-6]);
gg=[zeros(1,3);zeros(1,3);eye(3)];
gQg=gg*Q*gg';
sigma=[x(6:10,1) x(11:15,1) x(16:20,1) x(21:25,1) x(26:30,1)];
xdot2=A*sigma+sigma*A'+gQg;
xdot(6:10,1)=xdot2(:,1);
xdot(11:15,1)=xdot2(:,2);
xdot(16:20,1)=xdot2(:,3);
xdot(21:25,1)=xdot2(:,4);
xdot(26:30,1)=xdot2(:,5);
end

%generating sigma points
function [X,xw]=sigmap(mex1,covx)

    L=length(mex1);
    mex=reshape(mex1,[L,1]);
    ki=0;
    A=sqrtm((L+ki)*covx);
    
    X=zeros(L,2*L+1);
    X(:,1)=mex;

    xw(1)=ki/(L+ki);

    for ii=2:(L+1)
        X(:,ii)=mex+A(:,ii-1);
        xw(ii)=1/(2*(L+ki));
    end    
    for ii=(L+2):(2*L+1)
        X(:,ii)=mex-A(:,ii-L-1);
        xw(ii)=1/(2*(L+ki));
    end  
end

%output function    
function h=out(b)
global R0
h(1,1)=((b(1,1)-R0)^2+b(2,1)^2)^0.5;
h(2,1)=atan(b(2,1)/(b(1,1)-R0));
end
