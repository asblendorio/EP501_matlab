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
        for k = 1:100 % add z dimension 
            
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

%gradient of scalar function
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
    gradz(:,:,iz)=(f(:,:,iz+1)-f(:,:,iz-1))/2/dz;    %\partial/\partial y
end %for
grady(:,:,lz)=(f(:,:,lz)-f(:,:,lz-1))/dz;

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
    divz(:,:,iz)=(s(:,:,iz+1)-s(:,:,iz-1))/2/dz;    %\partial/\partial z
end %for
divz(:,:,lz)=(s(:,:,lz)-s(:,:,lz-1))/dz;

laplacian=divx+divy+divz;    %this is really laplacian b/c input is gradient

%% Integration using the Iterated Trapezoidal method
%constants
Q=1;
a=1;
e0=8.854.*10.^(-12);

%gradient of scalar function
dx=x(2)-x(1);
dy=y(2)-y(1);
dz=z(2)-z(1);
  
q = (e0.*laplacian).*f;
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

% Analytical Laplacian in 3D:
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

