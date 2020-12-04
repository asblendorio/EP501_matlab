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
dphi_dx_a = 1000; %volts

ep = e0*(10*tanh((x-dx)/l)-(10*tanh((x-dxx)/l)));
  
figure(1);
plot(ep,'b--');
xlabel('1/m');
ylabel('Farads');
title('Dielectric Function as a function of x');

%% 1B
%second order, centered
dy_dx=zeros(lx,1);
dx=x(2)-x(1);

%forward difference at the beginning
dy_dx(1)=(y(2)-y(1))/dx;

%centered difference on the interior
for ix=2:lx-1
    dy_dx(ix)=(y(ix+1)-y(ix-1))/2/dx;
end %for

%backward difference at the end
dy_dx(lx)=(y(lx)-y(lx-1))/dx;

figure(2);
plot(x,dy_dx,'k--');
xlabel('fill-in');
ylabel('fill-in');
title('fill-in');


%% 1C

% e = zeros(100,100);
% 
% for i=1:lx
%     e(x(i)) = e0.*(10.*tanh((x(i)-dx(i))./l)-(10.*tanh((x(i)-dxx(i))./l)));
% end %for  

