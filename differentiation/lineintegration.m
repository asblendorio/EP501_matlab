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
figure;
pcolor(x,y,r_phi);
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

% Plotting Magnetic Field B in terms of r 

for i=1:r_phi
    
end %for    
figure(2);
pcolor(x,y,r_phi);
shading flat;
colorbar;
hold off;



%% Analytical derivative 
dx = -phi_grid.*(r0).*sin(phi_grid);
dy = phi_grid.*(r0).*cos(phi_grid);
figure;
plot(dx,dy,'color','black','LineWidth',1.25)    
shading flat;
colorbar;
hold off;

%% Numerical Derivative

