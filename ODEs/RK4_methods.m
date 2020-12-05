%% Gridding in time
N=25;
tmin=0;
tmax=10;
t=linspace(tmin,tmax,N);
dt=t(2)-t(1);


%% Test problem, true solution
y0=1;
x0=1;
alpha=2;
ybar=y0*exp(-alpha*t);
xbar=x0*exp(-alpha*t);

%% Second order method; RK2
yRK2=zeros(1,N);
yRK2(1)=y0;
for n=2:N
    yhalf=yRK2(n-1)+dt/2*(-alpha*yRK2(n-1));
    yRK2(n)=yRK2(n-1)+dt*(-alpha*yhalf);
end %for

xRK2=zeros(1,N);
xRK2(1)=x0;
for n=2:N
    xhalf=xRK2(n-1)+dt/2*(-alpha*xRK2(n-1));
    xRK2(n)=xRK2(n-1)+dt*(-alpha*xhalf);
end %for

%% RK4 example; comparison against first and second order methods
yRK4=zeros(1,N);
yRK4(1)=y0;
for n=2:N
    dy1=dt*fRK(t(n-1),yRK4(n-1),alpha);
    dy2=dt*fRK(t(n-1)+dt/2,yRK4(n-1)+dy1/2,alpha);
    dy3=dt*fRK(t(n-1)+dt/2,yRK4(n-1)+dy2/2,alpha);
    dy4=dt*fRK(t(n-1)+dt,yRK4(n-1)+dy3,alpha);
    
    yRK4(n)=yRK4(n-1)+1/6*(dy1+2*dy2+2*dy3+dy4);
end %for

% x-component 
xRK4=zeros(1,N);
xRK4(1)=x0;
for n=2:N
    dx1=dt*fRK(t(n-1),xRK4(n-1),alpha);
    dx2=dt*fRK(t(n-1)+dt/2,xRK4(n-1)+dx1/2,alpha);
    dx3=dt*fRK(t(n-1)+dt/2,xRK4(n-1)+dx2/2,alpha);
    dx4=dt*fRK(t(n-1)+dt,xRK4(n-1)+dx3,alpha);
   
    xRK4(n)=xRK4(n-1)+1/6*(dx1+2*dx2+2*dx3+dx4);
end %for


%% Plots of RK solutions against true solution
figure(1);
clf;
plot(t,ybar,'o-');
xlabel('t');
ylabel('y(t)');
set(gca,'FontSize',20);
figure(1);
hold on;
plot(t,yRK2,'--')
figure(1);
plot(t,yRK4,'^-')
legend('exact','RK2','RK4')


%% RK2 stability considerations, FDE analysis
adt=linspace(0.01,3,20);
ladt=numel(adt);
G=zeros(ladt,1);
for igain=1:ladt
    G(igain)=(1-adt(igain)+1/2*adt(igain).^2);
end %for
figure(2);
plot(adt,G,'o')
set(gca,'FontSize',20);
xlabel('\alpha \Delta t');
ylabel('gain factor');