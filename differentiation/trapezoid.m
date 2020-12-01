%% Multidimensional function and partial derivatives:  grad and div (laplacian and curl also useful)
lx=100;
ly=100;
% lz=100;

%constants
Q=1;
a=1;
e0=8.854.*10.^(-12);
C = Q/(4*pi*e0*a);
phi=zeros(100,100);

x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);
% z = linspace(-3*a,3*a,lz);
[X,Y]=meshgrid(x,y);

dp_dx=zeros(size(phi));
dp_dy=zeros(size(phi));
% dp_dz=zeros(size(phi));

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
figure(1);
pcolor(x,y,f);
shading flat;

figure(2);
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
% gradx=zeros(size(f));
% grady=zeros(size(f));
% gradz=zeros(size(f));

%x component of gradient
dp_dx(:,1)=(phi(:,2)-f(:,1))/dx;
for ix=2:lx-1
    dp_dx(:,ix)=(phi(:,ix+1)-phi(:,ix-1))/2/dx;    %\partial/\partial x
end %for
dp_dx(:,lx)=(f(:,lx)-f(:,lx-1))/dx;

%y component of gradient
dp_dy(1,:)=(phi(2,:)-phi(1,:))/dy;
for iy=2:ly-1
    dp_dy(iy,:)=(phi(iy+1,:)-phi(iy-1,:))/2/dy;    %\partial/\partial y
end %for
dp_dy(ly,:)=(phi(ly,:)-phi(ly-1,:))/dy;

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

figure(4);
surface(x,y,div);
set(gca,'FontSize',24);
xlabel('x');
ylabel('y');
zlabel('z');
title('laplacian(f)');
colorbar;

%% Integration in multiple dimensions
%using an iterated trapezoidal method (sweeps of single dimensional integrations) or multi-dimensional
%trapezoidal method program that you write.


%terms inside integral
int=(e_0*laplace_phi).*phi;
s=zeros(size(Y,1));


%Integrating F at every y and z for all values of x 
for k=1:size(Z,1)-1
    for j=1:size(Y,3)-1 
        for i=2:size(X,2)-1
        s(j,k)=s(j,k)+(2*int(j,i,k));
        end
        s(j,k)=s(j,k)+int(j,1,k)+int(j,size(X,2),k);
    end
end

plane=zeros(size(Z,3),1);
%Integration for values of y at every z 
for k=1:size(Z,1)
    for j=2:size(Y,3)-1
    plane(k)=plane(k)+(2*surface(j,k));
    end
    plane(k)=plane(k)+surface(1,k)+surface(size(Y,3),k);             
    %Integration for last and first term
end

%Integrating for all values of z 
for k=2:size(Z,1)-1
    Integrated_term=Integrated_term+(2*plane(k));
end

Integrated_term=Integrated_term+plane(1)+plane(size(Z,3)); We=(-0.5)*Integrated_term;
fprintf('\n Total electrostatic energy We = %f Joules\n',We);
Total electrostatic energy We = 4978705645790294.000000 Joules
