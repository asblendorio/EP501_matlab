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
%f=exp(-(X.^2)/2/2).*exp(-Y.^2/2/1);

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
% z = linspace(-3*a,3*a,lz);
[X,Y]=meshgrid(x,y);
%f=exp(-(X.^2)/2/2).*exp(-Y.^2/2/1);

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

figure(6);
contourf(x,y,f);
xlabel('x');
ylabel('y');
title('f(x,y) and grad(f)');
colorbar;

%gradient of scalar function
dx=x(2)-x(1);
dy=y(2)-y(1);
% dz=z(2)-z(1);

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

%z component of gradient
% gradz(1,:)=(f(2,:)-f(1,:))/dz;
% for iz=2:lz-1
%     gradz(iz,:)=(f(iz+1,:)-f(iz-1,:))/2/dz;    %\partial/\partial y
% end %for
% gradz(lz,:)=(f(lz,:)-f(lz-1,:))/dz;

%add quiver on top of color plot
hold on;
quiver(X,Y,gradx,grady,'Color','white','LineWidth',2);
set(gca,'FontSize',24);
hold off;

%% Take the Laplacian by taking divergence of the previously computed gradient
f=gradx;
g=grady;
% h=gradz;

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


% z-derivative part of the divergence 
% divz=zeros(size(z));
% divz(1,:)=(h(2,:)-h(1,:))/dz;
% for iz=2:lz-1
%     divz(iz,:)=(h(iz+1,:)-h(iz-1,:))/2/dz;
% end %for
% divz(lz,:)=(h(lz,:)-h(lz-1,:))/dz;

div=divx+divy;    %this is really laplacian b/c input is gradient

figure(7);
surface(x,y,div);
set(gca,'FontSize',24);
xlabel('x');
ylabel('y');
zlabel('z');
title('laplacian(f)');
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
phi=zeros(100,100,100);

x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);
z = linspace(-3*a,3*a,lz);
[X,Y,Z]=meshgrid(x,y,z);

% dp_dx=zeros(size(phi));
% dp_dy=zeros(size(phi));
% dp_dz=zeros(size(phi));

%gradient of scalar function
dx=x(2)-x(1);
dy=y(2)-y(1);
dz=z(2)-z(1);
DE = 1; % energy initiliaztion 

for i=1:50
    for j=1:50
        for k=1:50
            if (sqrt(x(i).^2 +y(j).^2+z(k).^2) < a)
            phi(i,j,k) = (-1/2.*e0.*-6.*(Q/(8.*pi.*e0.*a.^3))).*(Q/(4.*pi.*e0.*a)-(Q/(8.*pi.*e0.*a.^3))*sqrt(x(i).^2 +y(j).^2+z(k).^2).^2-a.^2);
            
            else
            phi(i,j,k) = 0; %boundary conditions sets this value at zero outside of source field
            end %if
            % Integral 
            DE = ((0.5.*phi(i+1,j,k)+phi(i,j,k).*dx)+ (0.5.*phi(i,j+1,k)+phi(i,j,k).*dy)+(0.5.*phi(i,j,k+1)+phi(i,j,k).*dz) + DE);
        end %for  
    end %for
end %for
% We = WE*(2.15*10.^-3);
fprintf('\n Total electrostatic energy We = %e Joules\n',We);

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

disp('%%%%%%%%%%%%%%%%%%PROBLEM #4 ANSWER END%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%PROBLEM #4 ANSWER PLOTS%%%%%%%%%%%%%%%%%%');
%% END PROJECT #5