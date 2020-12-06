%% RK2 and systems of equations, oscillating charge example
% charged particle moving in the x-y direction in a magnetic field which
% has only a z-component...
%
% m dv/dt = q v x B 
% (ma = F)
%
% v = (vx,vy); B=(0,0,B);
%
% Resulting in the following system of equations:
%    m dvx/dt = q vy B
%    m dvy/dt = -q vx B

q=-1.6e-19;
m=1.67e-27;
B=5e-6;
omega=q*B/m;    %frequency of oscillation (can be shown via solution by hand gives a SHO)
tmin=0;
tmax=2*2*pi/abs(omega);    % follow particle for one oscillation periods
t=linspace(tmin,tmax,100);
dt=t(2)-t(1);
lt=numel(t);

% for j=75:-1:1
%     omega=q*B/m;    %frequency of oscillation (can be shown via solution by hand gives a SHO)
%     tmin(j)=0;
%     tmax(j)=2*2*pi/abs(omega);    % follow particle for one oscillation periods
%     t(j)=linspace(tmin,tmax,75);
%     dt(j)=t(2)-t(1);
%     lt=numel(t); 
%     if 
%         
%     end %    
% end %for    

% RK2 Method
vx=zeros(1,lt);
vy=zeros(1,lt);
vx(1)=1;     % vx initial conditions
vy(1)=1;     % vy initial conditions
% Loop for applying RK2 to a system of two equations
for n=2:lt
    %step x and y components together, this is the half update
    vxhalf=vx(n-1)+dt/2*(omega*vy(n-1));
    vyhalf=vy(n-1)-dt/2*(omega*vx(n-1));
    
    %now the full update
    vx(n)=vx(n-1)+dt*(omega*vyhalf);
    vy(n)=vy(n-1)-dt*(omega*vxhalf);    
end %for


%RK4 method
vx4=zeros(1,lt);
vy4=zeros(1,lt);
vx4(1)=1e3;     % vx initial conditions
vy4(1)=1e3;     % vy initial conditions

% Loop for applying RK4 to a system of two equations
for n=2:lt
    k1x=dt*(omega*vy4(n-1));    %k1 for the x differential equation
    k1y=-dt*(omega*vx4(n-1));    %k1 for the y differential equation
    
    k2x=dt*omega*(vy4(n-1)+k1y/2);
    k2y=-dt*omega*(vx4(n-1)+k1x/2);
    
    k3x=dt*omega*(vy4(n-1)+k2y/2);
    k3y=-dt*omega*(vx4(n-1)+k2x/2); 
    
    k4x=dt*omega*(vy4(n-1)+k3y);
    k4y=-dt*omega*(vx4(n-1)+k3x); 
    
    vx4(n)=vx4(n-1)+1/6*(k1x+2*k2x+2*k3x+k4x);
    vy4(n)=vy4(n-1)+1/6*(k1y+2*k2y+2*k3y+k4y);
end %for

%RK2 Integrate velocity to get position as a fn. of time, this assumes that the
%particle is initially at x,y = (0,0)
x=cumtrapz(t,vx);    %Matlab built-in for accumulating an integral value
y=cumtrapz(t,vy);
vz=1e3;
z=vz*t;

%RK4 integrate velocity to get position as a fn. of time, this assumes that the
%particle is initially at x,y = (0,0)
x4=cumtrapz(t,vx4);
y4=cumtrapz(t,vy4);
vz4=1e3;
z4=vz4*t;

% Plot velocity solutions for both RK2 and RK4
figure(1);
ax=plotyy(t,vx,t,vy);
set(ax(1),'FontSize',20);
set(ax(2),'FontSize',20);
xlabel('time (s)');
ylabel(ax(1),'v_x');
ylabel(ax(2),'v_y');
title('Runge-Kutta 2nd Order Method');

figure(2);
ax4=plotyy(t,vx4,t,vy4);
set(ax(1),'FontSize',20);
set(ax(2),'FontSize',20);
xlabel('time (s)');
ylabel(ax(1),'v_x');
ylabel(ax(2),'v_y');
title('Runge-Kutta 4th Order Method');

% Comet plot demo
figure(3);
comet3(x,y,z);
set(gca,'FontSize',20);
xlabel('x');
ylabel('y');
zlabel('z');
title('Runge-Kutta 2nd Order Method');

figure(4);
comet3(x4,y4,z4);
set(gca,'FontSize',20);
xlabel('x');
ylabel('y');
zlabel('z');
title('Runge-Kutta 4th Order Method');

figure(8);
ax=plotyy(t,vy,t,vy4);
set(ax(1),'FontSize',20);
set(ax(2),'FontSize',20);
xlabel('time (s)');
ylabel(ax(1),'v_x');
ylabel(ax(2),'v_y');
title('Comparison of Runge-Kurtta 2nd and 4th Order Methods');
