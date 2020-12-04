%% Author: Alec Sblendorio 
%% Due Date: December 4, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% Project #5
%%Problem 1: Numerical Vector Derivatives (curl)
%%Problem #2: Numerical Vector Derivatives (gradient and laplacian)
%%Problem #3: Integration in Multiple Dimensions
%%Problem #4: Line Integration
% Collaboration with Dennis Turano, and Bartosz Kwiecisnki, and Joseph
% Akash

%% Problem #1: Numerical Vector Derivatives (curl):
%%Part A: Plot the two components of the vector magnetic field defined by the 
%%piecewise function. Use an image plot (e.g. pcolor and shading flat in
% MATLAB for each magnetic field component (Bx,By) and have your plot 
%%show the region −3a ≤ x ≤ 3a, −3a ≤ y ≤ 3a. Make sure you add a colorbar 
%%and axis labels to your plot. You will need to define a range and resolution
%%in x and y, and create a meshgrid from that. Be sure to use a resolution 
%%fine enough to resolve important variations in this function.

%%Part B: Make a quiver plot of the magnetic field B; add labels, etc

%%Part C: Compute the numerical curl of B, i.e. ∇ × B. Use centered 
%%differences on the interior grid points and first-order derivatives on the edges.
%%Plot your result using imagesc, or pcolor.

%%Part D: Compute ∇ × B analytically (viz. by hand). Plot the alongside your
%%numerical approximation and demonstrate that they are suitably similar.

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');

%% Multidimensional function and partial derivatives:  grad and div (laplacian and curl also useful)
lx=100;
ly=100;
%constants
I=10;
a=0.005;
m0=(4*pi)*10.^(-7);

x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);
[X,Y]=meshgrid(x,y);

f = zeros(100,100);
for i=1:100
    for j=1:100
        if (sqrt(x(i).^2 +y(j).^2) < a)
            f(i,j) = ((m0.*I)/(2.*pi.*a.^2).*sqrt(x(i).^2+y(j).^2)).*(-y(j)/sqrt(x(i).^2+y(j).^2)+(x(i)./sqrt(x(i).^2+y(j).^2)));
        else
            f(i,j) = ((m0.*I)/(2.*pi.*sqrt(x(i).^2+y(j).^2))).*(-y(j)/sqrt(x(i).^2+y(j).^2)+(x(i)./sqrt(x(i).^2+y(j).^2)));
        end %if 
    end %for
end %for

%% plotter
figure(1);
pcolor(x,y,f);
shading flat;
colorbar;
hold off;

figure(2);
contourf(x,y,f);
xlabel('x');
ylabel('y');
xlim([-0.015 0.015]);
ylim([-0.015 0.015]);
title('f(x,y) and grad(f)');
colorbar;
hold off;

%gradient of scalar function
dx=x(2)-x(1);
dy=y(2)-y(1);

%% Gradient
gradx=zeros(size(f));
grady=zeros(size(f));

%x component of gradient
gradx(:,1)=(f(:,2)-f(:,1))/dx;
for ix=2:lx-1
    gradx(:,ix)=(f(:,ix+1)-f(:,ix-1))/2/dx;    %\partial/\partial x
end %for
gradx(:,lx)=(f(:,lx)-f(:,lx-1))/dx;

%y component of gradient
grady(1,:)=(f(2,:)-f(1,:))/dy;
for iy=2:ly-1
    grady(iy,:)=(f(iy+1,:)-f(iy-1,:))/2/dy;    %\partial/\partial y
end %for
grady(ly,:)=(f(ly,:)-f(ly-1,:))/dy;

%add quiver on top of color plot
hold on;
quiver(X,Y,gradx,grady,'Color','white','LineWidth',2);
set(gca,'FontSize',24);


%% Curl Numerical
curlx=zeros(ly,lx);
curly=zeros(ly,lx);

% x component of curl
curlx(:,1)=(f(:,2)-f(:,1))/dx;
for ix=2:lx-1
    curlx(:,ix)=(f(:,ix+1)-f(:,ix-1))/2/dx;    %\partial/\partial x
end %for
curlx(:,lx-1)=(f(:,lx-1)-f(:,lx-2))/dx;

%y component of curl 
curly(1,:)=(f(2,:)-f(1,:))/dy;
for iy=2:ly-1
    curly(iy,:)=(f(iy+1,:)-f(iy-1,:))/2/dy;    %\partial/\partial y
end %for
curly(ly-1,:)=(f(ly-1,:)-f(ly-2,:))/dy;

curl_B=curlx-curly;

figure(3);
imagesc(x,y,curl_B);
xlabel('x');
ylabel('y');
title('Curl of B');
colorbar;
hold off;

%% Curl Analytical
f = zeros(100,100);
for i=1:100
    for j=1:100
        if (sqrt(x(i).^2 +y(j).^2) < a)
            f(i,j) = 2*((m0.*I)/(2.*pi.*a.^2)); %Hand calculated curl 
        else
            f(i,j) = 0; %Hand calculated curl 
        end %if 
    end %for
end %for

figure(4);
imagesc(x,y,f);
xlabel('x');
ylabel('y');
title('Analytical Curl of B');
colorbar;
shading flat; 


disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER END%%%%%%%%%%%%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER PLOTS%%%%%%%%%%%%%%%%%%');
%% Problem #2: Numerical Vector Derivatives (gradient and laplacian)
%%Part A:  Compute and plot the scalar field. and plot this function in the 
%%region−3a≤x≤3a,−3a≤y≤3a in the z=0 plane. Be sure to use a resolution 
%%fine enough to resolve variations in this function 
%%(aside from those associated with the singularity).

%%Part B: Write a function to numerically compute the Laplacian of a scalar field, i.e. ∇2Φ.
%%Plot your result with appropriate labels and colorbars.

%%Part C: Compute an analytical laplacian (viz. differentiate by hand),
%%plot the results alongside your numerical calculation, 
%%and demonstrate that your numerical laplacian is suitably accurate.

disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
%% Multidimensional function and partial derivatives:  grad and div (laplacian and curl also useful)
lx=100;
ly=100;
% lz=100;

%constants
Q=1;
a=1;
e0=8.854.*10.^(-12);

x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);
[X,Y]=meshgrid(x,y);

f = zeros(100,100);
for i=1:100
    for j=1:100
        if (sqrt(x(i).^2 +y(j).^2) < a)
            f(i,j) = Q/(4.*pi.*e0.*a) - Q/(8.*pi.*e0.*a.^3).*(x(i).^2+y(j).^2-a.^2); 
        else
            f(i,j) = Q/(4.*pi.*e0.*sqrt(x(i).^2 +y(j).^2));
        end %if 
    end %for
end %for

%% plotter
figure(5);
pcolor(x,y,f);
xlabel('x');
ylabel('y');
title('Scalar Field');
shading flat;
colorbar;

figure(6);
contourf(x,y,f);
xlabel('x');
ylabel('y');
title('f(x,y) and grad(f)');
colorbar;

%gradient of scalar function
dx=x(2)-x(1);
dy=y(2)-y(1);

%% Gradient
gradx=zeros(size(f));
grady=zeros(size(f));
% gradz=zeros(size(f));

%x component of gradient
gradx(:,1)=(f(:,2)-f(:,1))/dx;
for ix=2:lx-1
    gradx(:,ix)=(f(:,ix+1)-f(:,ix-1))/2/dx;    %\partial/\partial x
end %for
gradx(:,lx)=(f(:,lx)-f(:,lx-1))/dx;

%y component of gradient
grady(1,:)=(f(2,:)-f(1,:))/dy;
for iy=2:ly-1
    grady(iy,:)=(f(iy+1,:)-f(iy-1,:))/2/dy;    %\partial/\partial y
end %for
grady(ly,:)=(f(ly,:)-f(ly-1,:))/dy;


%add quiver on top of color plot
hold on;
quiver(X,Y,gradx,grady,'Color','white','LineWidth',2);
set(gca,'FontSize',24);
hold off;

%% Take the Laplacian by taking divergence of the previously computed gradient
f=gradx;
g=grady;

%x-derivative part of the divergence
divx=zeros(size(f));
divx(:,1)=(f(:,2)-f(:,1))/dx;
for ix=2:lx-1
    divx(:,ix)=(f(:,ix+1)-f(:,ix-1))/2/dx;
end %for
divx(:,lx)=(f(:,lx)-f(:,lx-1))/dx;

%y-derivative part of the divergence
divy=zeros(size(y));
divy(1,:)=(g(2,:)-g(1,:))/dy;
for iy=2:ly-1
    divy(iy,:)=(g(iy+1,:)-g(iy-1,:))/2/dy;
end %for
divy(ly,:)=(g(ly,:)-g(ly-1,:))/dy;

laplacian=divx+divy;    %this is really laplacian b/c input is gradient

figure(7);
surface(x,y,laplacian);
set(gca,'FontSize',24);
xlabel('x');
ylabel('y');
zlabel('z');
title('Numerically Computed Laplacian(f)');
colorbar;
shading flat; 
disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER END%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER PLOTS%%%%%%%%%%%%%%%%%%');
%% Problem #3: Integration in Multiple Dimensions.
%%Part A:  Numerically compute the electrostatic energy in the region
%%R ≡ −3a ≤ x ≤ 3a,−3a ≤ y ≤ 3a, −3a ≤ z ≤ 3a, defined by the integral.
%%using an iterated trapezoidal method (sweeps of single dimensional integrations) 
%%or multi-dimensional trapezoidal method program that you write.

disp('%%%%%%%%%%%%%%%%%%PROBLEM #3 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
%% Integration using the Iterated Trapezoidal method
lx=100;
ly=100;
lz=100;

%constants
Q=1;
a=1;
e0=8.854.*10.^(-12);

x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);
z = linspace(-3*a,3*a,lz);
[X,Y,Z]=meshgrid(x,y,z);

f = zeros(100,100,100);

for i=1:100
    for j=1:100
        for k = 1:100 % add 3rd dimension 
            
        if (sqrt(x(i).^2 +y(j).^2+z(k).^2) < a)
            f(i,j,k) = Q/(4.*pi.*e0.*a) - Q/(8.*pi.*e0.*a.^3).*(x(i).^2+y(j).^2+z(k).^2-a^2); 
        else
            f(i,j,k) = Q/(4.*pi.*e0.*sqrt(x(i).^2 +y(j).^2+z(k).^2));
        end %if
        
        end %for
    end %for
end %for

% dp_dx=zeros(size(phi));
% dp_dy=zeros(size(phi));
% dp_dz=zeros(size(phi));

%% Gradient of scalar function
dx=x(2)-x(1);
dy=y(2)-y(1);
dz=z(2)-z(1);

gradx=zeros(size(f));
grady=zeros(size(f));
gradz=zeros(size(f));

% x component of gradient
gradx(:,1,:)=(f(:,2,:)-f(:,1,:))/dx;
for ix=2:lx-1
    gradx(:,ix,:)=(f(:,ix+1,:)-f(:,ix-1,:))/2/dx;    %\partial/\partial x
end %for
gradx(:,lx,:)=(f(:,:,lx)-f(:,:,lx-1))/dx;

% y component of gradient
grady(1,:,:)=(f(2,:,:)-f(1,:,:))/dy;
for iy=2:ly-1
    grady(iy,:,:)=(f(iy+1,:,:)-f(iy-1,:,:))/2/dy;    %\partial/\partial y
end %for
grady(ly,:,:)=(f(ly,:,:)-f(ly-1,:,:))/dy;

% z component of gradient
gradz(:,:,1)=(f(:,:,2)-f(:,:,1))/dz;
for iz=2:lz-1
    gradz(:,:,iz)=(f(:,:,iz+1)-f(:,:,iz-1))/2/dz;    %\partial/\partial z
end %for
grady(:,:,lz)=(f(:,:,lz)-f(:,:,lz-1))/dz;

%% Divergence of Gradient
k=gradx;
g=grady;
s=gradz;
 
% x-derivative part of the divergence
divx=zeros(size(k));
divx(:,1,:)=(k(:,2,:)-k(:,1,:))/dx;
for ix=2:lx-1
    divx(:,ix,:)=(k(:,ix+1,:)-k(:,ix-1,:))/2/dx;
end %for
divx(:,lx,:)=(k(:,lx,:)-k(:,lx-1,:))/dx;

% y-derivative part of the divergence
divy=zeros(size(g));
divy(1,:,:)=(g(2,:,:)-g(1,:,:))/dy;
for iy=2:ly-1
    divy(iy,:,:)=(g(iy+1,:,:)-g(iy-1,:,:))/2/dy;
end %for
divy(ly,:,:)=(g(ly,:,:)-g(ly-1,:,:))/dy;

% z-derivative part of the divergence 
divz=zeros(size(s));
divz(:,:,1)=(s(:,:,2)-s(:,:,2))/dz;
for iz=2:lz-1
    divz(:,:,iz)=(s(:,:,iz+1)-s(:,:,iz-1))/2/dz;    
end %for
divz(:,:,lz)=(s(:,:,lz)-s(:,:,lz-1))/dz;

laplacian2=divx+divy+divz;    %this is really laplacian b/c input is gradient

%% Integration using the Iterated Trapezoidal method
%gradient of scalar function
dx=x(2)-x(1);
dy=y(2)-y(1);
dz=z(2)-z(1);
% input numerical laplacian  
q = (e0.*laplacian2).*f;
DE = 1; % energy initiliaztion 

for i=1:lx-1
    for j=1:ly-1
        for k=1:lz-1
            
            % Integral 
            integral = (-0.0625)*...
            (q(i,j,k) + q(i+1,j,k)+ ...
            (q(i,j+1,k) + q(i+1,j+1,k)+ ...
            (q(i,j,k+1) + q(i+1,j,k+1)+...
            (q(i,j+1,k+1) + q(i+1,j+1,k+1)))));
            DE = DE + integral.*dx.*dy.*dz;
            
        end %for  
    end %for
end %for

We = DE;
fprintf('\n Total electrostatic energy We (Numerical) = %e Joules\n',We);

%% Analytical Laplacian:
an_lap = zeros(100,100,100);


for i = 1:lx
    for j = 1:ly
        for k = 1:lz
            
        if sqrt(x(i).^2+y(j).^2+z(k).^2) < a
          an_lap(i,j,k) = -6*(Q/(8.*pi.*e0.*a.^3)); 
        else 
          an_lap(i,j,k) = 0;
        end %if
        
        end %for
    end %for
end %for

q2 = (e0.*an_lap).*f;
DE2 = 1; % energy initiliaztion 

for i=1:lx-1
    for j=1:ly-1
        for k=1:lz-1
            % Integral 
            integral2 = (-0.0625)*...
            (q2(i,j,k) + q2(i+1,j,k)+ ...
            (q2(i,j+1,k) + q2(i+1,j+1,k)+ ...
            (q2(i,j,k+1) + q2(i+1,j,k+1)+...
            (q2(i,j+1,k+1) + q2(i+1,j+1,k+1)))));
           DE2 = DE2 + integral2.*dx.*dy.*dz;
        end %for  
    end %for
end %for

We2 = DE2;
fprintf('\n Total electrostatic energy We (Analytical) = %e Joules\n',We2);

disp('%%%%%%%%%%%%%%%%%%PROBLEM #3 ANSWER END%%%%%%%%%%%%%%%%%%');
%% Problem #4: Line Integration
%%Part A:  Compute and plot the parametric path
%%r(φ)≡x(φ)eˆ +y(φ)eˆ =r cosφeˆ +r sinφeˆ (0≤φ≤2π) xy0x0y
%%in the x,y plane on the same axis as your magnetic field components 
%%from problem 1 (create a new figure which plots the path on top of a 
%%pcolor plot of each component). Take r0 = 2a. You will need to define a 
%%grid in φ to do this.

%%Part B: Plot the two components of the magnetic field B(x(φ), y(φ)) at 
%%the x, y points along r and visually compare against your image plots of
%%the magnetic field and path to verify.

%%Part C: Numerically compute the tangent vector to the path r by performing the derivative:
%%Compare your numerical results against the analytical derivative (e.g. plot the two)
%%and refine your grid in φ (if necessary) such that you get visually 
%%acceptable results - i.e. such that the path appears circular.

%%Part D: Numerically compute the auxiliary magnetic field integrated around the path r, i.e.
%%where the differential path length is given by:

disp('%%%%%%%%%%%%%%%%%%PROBLEM #4 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
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
figure(8);
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
figure(9);
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
figure(10);
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
    B_phi1(i) = fprime1.*(Bx_phi(i)+By_phi1(i));
    
end %for 

integral3 = 0;
for j=1:99
    integral3 = (0.25).*((B_phi1(i+1)+B_phi1(i).*dx))+integral3;
end %for    

fprintf('\n The Magnetic field has a current of = %f Amps\n',integral3);

disp('%%%%%%%%%%%%%%%%%%PROBLEM #4 ANSWER END%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%PROBLEM #4 ANSWER PLOTS%%%%%%%%%%%%%%%%%%');
%% END PROJECT #5