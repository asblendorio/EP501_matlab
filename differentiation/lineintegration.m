%% Line Integration 
% Compute and plot the parametric path 
lx=100;
ly=100;
%constants
I=10;
a=0.005;
m0=(4*pi)*10.^(-7);
r0=2*a;

x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);
[X,Y]=meshgrid(x,y);

r_phi = zeros(100,100);

%Magnetic Field component from Problem #1 
for i=1:100
    for j=1:100
        if (sqrt(x(i).^2 +y(j).^2) < a)
            r_phi(i,j) = ((m0.*I)/(2.*pi.*a.^2).*sqrt(x(i).^2+y(j).^2)).*(-y(j)/sqrt(x(i).^2+y(j).^2)+(x(i)./sqrt(x(i).^2+y(j).^2)));
        else
            r_phi(i,j) = ((m0.*I)/(2.*pi.*sqrt(x(i).^2+y(j).^2))).*(-y(j)/sqrt(x(i).^2+y(j).^2)+(x(i)./sqrt(x(i).^2+y(j).^2)));
        end %if 
    end %for
end %for
figure(1);
pcolor(x,y,r_phi);
xlabel('x');
ylabel('y');
title('Parametric Path vs. Magnetic Field ');
shading flat;
colorbar;
hold on;

% r_phix = zeros(1,100);
% r_phiy = zeros(1,100);

phi_grid = linspace(0,2*pi,lx);
x = r0.*cos(phi_grid);
y = r0.*sin(phi_grid);
%% Plotter
plot(x,y,'color','black','LineWidth',1.25)    
shading flat;
colorbar;
hold off;

%% Plotting Magnetic Field B in terms of r 
B_phi = zeros(1,100); 
Bx_phi= zeros(1,100); 
By_phi= zeros(1,100);

for i=1:lx
        B_phi(i) = ((m0*I)/(2*pi*(sqrt(x(i)^2 + ...
        y(i)^2)))*(-y(i)/(sqrt(x(i)^2 + y(i)^2)) + ...
        x(i)/(sqrt(x(i)^2 + y(i)^2))));
        
        Bx_phi(i)=((m0*I)/(2*pi*(sqrt(x(i).^2+y(i).^2))* ...
        (-y(i)./sqrt(x(i).^2+y(i).^2))));
        
        By_phi(i)=((m0*I)/(2*pi*(sqrt(x(i).^2+y(i).^2))* ...
        (x(i)./sqrt(x(i).^2+y(i).^2)))); 
    
end %for
figure(2);
plot(phi_grid,B_phi,'r','LineWidth',1);
xlabel('\Phi (in radians)'); 
ylabel('B_{x}(x,y)=B_{x}(\Phi) (in Tesla)'); 
title('Magnetic Field Components at r=0.01 m '); 
legend('B');
grid on;
hold off;


%% Numerical Derivative
phi_grid = linspace(0,2*pi,lx);
x = r0.*cos(phi_grid);
y = r0.*sin(phi_grid);

%second order, centered
dy_dx=zeros(lx,1);
dx=x(2)-x(1);
dx_dy=zeros(ly,1);
dy=y(2)-y(1);

%centered difference on the interior
for ix=2:lx-1
    dy_dx(ix)=(x(ix)-x(ix-1))/dx;
end %for
dy_dx(1)=dy_dx(2);
%backward difference at the end

%forward difference at the beginning
for iy=2:ly-1
    dx_dy(iy)=(y(iy)-y(iy-1))/dy;
end %for
dx_dy(1)=dx_dy(2);


%% Analytical derivative 
an_dx = -(r0).*sin(phi_grid);
an_dy = (r0).*cos(phi_grid);
%% Plotter
figure(3);
plot(an_dx,an_dy,'color','red','LineWidth',1.25);
xlabel('dx');
ylabel('f(x)');
title('Analytical Vs. Numerical Tangent Vector');
grid on;
hold on;
xlim([-3*a 3*a]);
ylim([-3*a 3*a]);

%% Auxiliary Magnetic Field 
B_phi1 = zeros(1,100); 
Bx_phi1= zeros(1,100); 
By_phi1= zeros(1,100);

for i=1:lx-1
    
    fprime1 = ((dy_dx(i)+dx_dy(i))./m0);
    By_phi1(i) = (r0).*cos(x(i));
    Bx_phi1(i) = -(r0).*sin(x(i));
    B_phi1(i) = fprime1.*(Bx_phi(i)+ By_phi1(i));
    
end %for 

integral3 = 0;
for j=1:99
    integral3 = (1/4).*((B_phi1(i+1)+B_phi1(i).*dx))+integral3;
end %for    

disp(round(integral3));


% dr_dphi = dx_dy+dy_dx;
% dphi = phi_grid(2)-phi_grid(1);
% % differential path length
% dl = dr_dphi*dphi;
% dl = sum(dl,'all');
% 
% B = sum(B_phi,'all');
% j = dl.*(B/m0);
% 
% % integral definition 
% current = (B/m0).*dl;
% current = sum(current,'all');
% xsoln = current./2;
% disp('The numerically computed auxiliary magnetic field integrated around the path r is:');
% disp(round(xsoln));


