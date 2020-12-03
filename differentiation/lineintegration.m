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

% phi_grid = linspace(0,2*pi,lx);
phi=(0:pi/100:2*pi);
x = r0.*cos(phi);
y = r0.*sin(phi);
%% Plotter
plot(x,y,'color','black','LineWidth',1.25)    
shading flat;
colorbar;
hold off;

% Plotting Magnetic Field B in terms of r 
rphi_x=r0*cos(phi); 
rphi_y=r0*sin(phi); 

rad=sqrt(rphi_x.^2+rphi_y.^2);

B_x_phi=zeros(size(x,1),1); 
B_y_phi=zeros(size(y,1),1);

for i=1:size(x,1)
    B_x_phi(i)=((-m0*I)/(2*pi*rad(i)*rad(i)))*y(i);
    B_y_phi(i)=((-m0*I)/(2*pi*rad(i)*rad(i)))*x(i);
end %for

B_phi=sqrt(B_x_phi.^2+B_y_phi.^2);

figure(2);
title('Plot of Magnetic field components along the pathline');
yyaxis left;
xlabel('\Phi (in radians)'); 
ylabel('B_{x}(x,y)=B_{x}(\Phi) (in Tesla)'); 
plot(phi,B_x_phi,'k-','LineWidth',2);
yyaxis right; 
plot(phi,B_y_phi,'r-','LineWidth',2);
hold on;

plot(phi,B_phi,'b--','LineWidth',1);
hold on;

legend('B_x','B_y','B (right axis)');
xlabel('\Phi (in radians)');
ylabel('B_{y}(x,y)=B_{y}(\Phi) (in Tesla)'); 
title('Magnetic field components plotted at r=0.01 m '); 
grid on;


%% Analytical derivative 
dx = -phi.*(r0).*sin(phi);
dy = phi.*(r0).*cos(phi);
% figure;
% plot(dx,dy,'color','black','LineWidth',1.25)    
% shading flat;
% colorbar;
% hold off;

%% Numerical Derivative

