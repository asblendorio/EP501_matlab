%% Author: Alec Sblendorio 
%% Due Date: December 6, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% Project #6
%%Problem 1: Electrostatic Potential in Dielectrics with Boundary
%%Conditions
%%Problem #2: Application of Ordinary Differential Equations
% Collaboration with Dennis Turano, and Bartosz Kwiecisnki, and Joseph
% Akash

%% Problem #1: Electrostatic Potential in Dielectrics with Boundary Conditions
%%Part A: Plot the dielectric function and note that it varies rapidly near the edges of the domain of interest.
%%Part B: Develop a system of finite difference equations for this system 
%%based on second order accurate centered differences (numerical differences
%%may also be used for dielectric function derivatives). 
%%Include this system in your homework submission.

%%Part C: Develop two finite difference equations for the boundary conditions
%%of this system. Use a first- order forward difference at the x = −a boundary.
%%Include these equations with your submission.

%%Part D: Solve your system of equations using the MATLAB “\” operator.
%%Part E: Since the dielectric function varies rapidly at the boundary, 
%%this is a problem where a second order (forward) difference may be useful
%%(see course notes for formula). Reformulate your matrix system to include
%%this for the x = −a boundary and solve the system numerically for 20 
%%grids points. Plot the result and compare it against the solution with a
%%first order forward difference.

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
%% Electrostatic potential of a dielectric function 
lx = 100;
ly = 100;
%constants 
e0=8.854.*10.^(-12);
a=0.01;
l=a./5;
x = linspace(-a,a,lx);
y = linspace(-a,a,ly);
[X,Y]=meshgrid(x,y);

dx = -9.*a./10; 
dxx = 9.*a./10;

%boundary conditions 
phi_a = 100; %volts
dphi_dx_a = 1000; %volts/m
bdcond = [dphi_dx_a;zeros(lx-2,1);phi_a]; 

ep = e0*(10*tanh((x-dx)/l)-(10*tanh((x-dxx)/l)));
  
figure(1);
plot(ep,'b--');
xlabel('1/m');
ylabel('Farads');
title('Dielectric Function as a function of x');

%% 1B
%second order, centered
dep_dx=zeros(lx,1);
dx=x(2)-x(1);
%forward difference at the beginning
dep_dx(1)=(ep(2)-ep(1))/2*dx;

%centered difference on the interior
for ix=2:lx-1
    dep_dx(ix)=(ep(ix+1)-ep(ix-1))/2/dx;
end %for

%backward difference at the end
dep_dx(lx)=(ep(lx)-ep(lx-1))/2*dx;

%%
dep_dx2=zeros(lx,1);
dxx=x(2)-x(1);
%first order derivative approximation (backward)
%interior
for ix=2:lx-1
    dep_dx2(ix)=(ep(ix)-ep(ix-1))/dxx;
end %for
dep_dx2(1)=dep_dx2(2);

%% Plotter
figure(2);
plot(x,dep_dx,'k--');
hold on;
plot(x,dep_dx2,'r--');
xlabel('x (m)');
ylabel('\epsilon(x)');
title('1st and 2nd order Finite Difference Derivative with %d grid points');
legend('2nd Order','1st Order');

%% 1D solve using Matlab Operator 
% ep2(:)=ep(:);
ep3=zeros(lx,ly);
ep3(1,:)=[-1/dx^2,1/dx,zeros(1,lx-2)];
ep3(lx,:)=[zeros(1,lx-1),1];

for i=2:lx-1
    ep3(i,i)=-2.*ep(i)./dx.^2;
    ep3(i,i+1)=ep(i)./dx.^2+dep_dx(i)./(2.*dx);
    ep3(i,i-1)=ep(i)./dx.^2-dep_dx(i)./(2.*dx);
end %for  

disp('Matlab,GNU/Octave built-in solution:');
xsoln=ep3\bdcond;
disp(xsoln);

figure(3);
plot(xsoln,'r--');
xlabel('x (m)');
ylabel('y');
title('Dielectric Function with x = -a at boundary');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER END%%%%%%%%%%%%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER PLOTS%%%%%%%%%%%%%%%%%%');
%% Problem #2: Application of Ordinary Differential Equations
%%Part A:  Solve this system with the RK4 algorithm and plot the velocity
%%and position as a function of time. Compare this against the RK2 solution
%%in the repository and show they are roughly equal. Use 2πm 100 time steps
%%per particle oscillation period T = qB .

%%Part B: Show that the RK4 solution is better in the sense that it can
%%solve the problem accurately with fewer time steps.

%%Part C: Suppose the magnetic field is changed to vary linearly in the y-direction:
%%B(y)=B 1+(1/2)y in e_z
%%Use your RK4 solver to find the velocity. Plot the velocities and path of 
%%the particle in the x-y plane for a least five periods of oscillation. 
%%HINT: The particle should execute trochoidal motion.

disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');

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
figure(4);
ax=plotyy(t,vx,t,vy);
set(ax(1),'FontSize',20);
set(ax(2),'FontSize',20);
xlabel('time (s)');
ylabel(ax(1),'v_x');
ylabel(ax(2),'v_y');
title('Runge-Kurtta 2nd Order Method');

figure(5);
ax4=plotyy(t,vx4,t,vy4);
set(ax(1),'FontSize',20);
set(ax(2),'FontSize',20);
xlabel('time (s)');
ylabel(ax(1),'v_x');
ylabel(ax(2),'v_y');
title('Runge-Kutta 4th Order Method');

% Comet plot demo
figure(6);
comet3(x,y,z);
set(gca,'FontSize',20);
xlabel('x');
ylabel('y');
zlabel('z');
title('Runge-Kutta 2nd Order Method');

figure(7);
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
title('Comparison of Runge-Kutta 2nd and 4th Order Methods');


%% Part C

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
figure(9);
ax4=plotyy(t,vx4,t,vy4);
xlabel('time (s)');
ylabel(ax(1),'v_x');
ylabel(ax(2),'v_y');
title('Runge-Kutta 4th Order Method');

% Comet plot demo
figure(10);
comet(x4,y4);
hold on;
set(gca,'FontSize',20);
xlabel('x(m)');
ylabel('y(m)');
legend('X(t)');
title('Runge-Kutta 4th Order Method');
disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER END%%%%%%%%%%%%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER PLOTS%%%%%%%%%%%%%%%%%%');


%% END PROJECT #6