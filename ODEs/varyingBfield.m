%% RK2 and RK4 systems of equations, oscillating charge example
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

q=1.6e-19;
m=1.67e-27;
B=5e-5;
omega=q*B/m;    %frequency of oscillation (can be shown via solution by hand gives a SHO)
tmin=0;
tmax=5*2*pi/abs(omega);    % follow particle for one oscillation periods
t=linspace(tmin,tmax,100);
dt=t(2)-t(1);
lt=numel(t);

%RK4 method
vx4=zeros(1,lt);
vy4=zeros(1,lt);
vx4(1)=1e3;     % vx initial conditions
vy4(1)=1e3;     % vy initial conditions
ynew=zeros(1,lt);
% Loop for applying RK4 to a system of two equations
for l=2:lt
    % set new value for the Magnetic Field each time as it varies in Y
    ynew(l) = ynew(l-1)+(dt*vy4(l-1));
    B2 = B*(1+.5*(ynew(l)));
    omega2 = q*B2/m;
    
    k1x=dt*(omega2*vy4(l-1));    %k1 for the x differential equation
    k1y=-dt*(omega2*vx4(l-1));    %k1 for the y differential equation
    
    k2x=dt*omega2*(vy4(l-1)+k1y/2);
    k2y=-dt*omega2*(vx4(l-1)+k1x/2);
    
    k3x=dt*omega2*(vy4(l-1)+k2y/2);
    k3y=-dt*omega2*(vx4(l-1)+k2x/2); 
    
    k4x=dt*omega2*(vy4(l-1)+k3y);
    k4y=-dt*omega2*(vx4(l-1)+k3x); 
    
    vx4(l)=vx4(l-1)+1/6*(k1x+2*k2x+2*k3x+k4x);
    vy4(l)=vy4(l-1)+1/6*(k1y+2*k2y+2*k3y+k4y);
  
end %for

%RK4 integrate velocity to get position as a fn. of time, this assumes that the
%particle is initially at x,y = (0,0)
x4=cumtrapz(t,vx4);
y4=cumtrapz(t,vy4);
vz4=1e3;
z4=vz4*t;

% Plot velocity solutions for  RK4
figure(2);
ax4=plotyy(t,vx4,t,vy4);
xlabel('time (s)');
ylabel('v_y');
title('Runge-Kutta 4th Order Method');

% Comet plot demo
figure(4);
comet(x4,y4);
hold on;
set(gca,'FontSize',20);
xlabel('x(m)');
ylabel('y(m)');
legend('X(t)');
title('Runge-Kutta 4th Order Method');