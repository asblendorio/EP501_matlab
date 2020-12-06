%% Electrostatic potential of a dielectric function 
lx = 20;
ly = 20;
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
title('1st and 2nd order Finite Difference Derivative');
legend('2nd Order','1st Order');

%% 1D solve using Matlab Operator 
ep2(:)=ep(:);
ep3=zeros(100,100);
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
