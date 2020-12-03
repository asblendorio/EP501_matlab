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
        
        Bx_phi(i)=((m0*I)/(2*pi*(sqrt(x(i).^2+y(i).^2))*...
        (-y(i)./sqrt(x(i).^2+y(i).^2))));
        
        By_phi(i)=((m0*I)/(2*pi*(sqrt(x(i).^2+y(i).^2))*...
        (x(i)./sqrt(x(i).^2+y(i).^2)))); 
    
end %for

figure(2);
title('Plot of Magnetic Field');
xlabel('\Phi (in radians)'); 
ylabel('B_{x}(x,y)=B_{x}(\Phi) (in Tesla)'); 

plot(phi_grid,B_phi,'b--','LineWidth',0.5);
hold on;
plot(phi_grid,Bx_phi,'r-','LineWidth',1);
hold on;
plot(phi_grid,By_phi,'k-','LineWidth',1);
hold on;

legend('Bx','By','B');
xlabel('\Phi'); 
ylabel('B_{y}(x,y)=B_{y}(\Phi) (in Tesla)'); 
title('Magnetic Field Components at r=0.01 m '); 
grid on;
hold off;

%% Analytical derivative 
dx = -(r0).*sin(phi_grid);
dy = (r0).*cos(phi_grid);
figure(3);
plot(dx,dy,'color','black','LineWidth',1.25);
xlabel('x');
ylabel('y');
title('Analytical vs. Numerical Tangent Vector');
shading flat;
colorbar;
hold off;

%% Numerical Derivative
%gradient of scalar function
dx=x(2)-x(1);
dy=y(2)-y(1);

dx_dphi=zeros(lx,1);


