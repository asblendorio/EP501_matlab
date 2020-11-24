%% Multidimensional function and partial derivatives:  grad and div (laplacian and curl also useful)
lx=100;
ly=100;
% lz=100;
%constants
I=10;
a=0.005;
m0=(4*pi)*10.^(-7);

% x=linspace(-5,5,lx);
% y=linspace(-5,5,ly);
x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);
% z = linspace(-3*a,3*a,lz);
[X,Y]=meshgrid(x,y);
% [X,Y,Z]=meshgrid(x,y,z);
 
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
hold off;

% %% Curl 
% curlx=zeros(size(f));
% curly=zeros(size(f));
% curlz=zeros(size(f));
% 
% %x component of curl
% curlx(:,1)=(f(:,2)-f(:,1))/dx;
% for ix=2:lx-1
%     curlx(:,ix)=(f(:,ix+1)-f(:,ix-1))/2/dx;    %\partial/\partial x
% end %for
% curlx(:,lx)=(f(:,lx)-f(:,lx-1))/dx;
% 
% %y component of curl 
% curly(1,:)=(f(2,:)-f(1,:))/dy;
% for iy=2:ly-1
%     curly(iy,:)=(f(iy+1,:)-f(iy-1,:))/2/dy;    %\partial/\partial y
% end %for
% curly(ly,:)=(f(ly,:)-f(ly-1,:))/dy;
% 
% %z component of curl 
% curlz(1,:)=(f(2,:)-f(1,:))/dz;
% for iz=2:lz-1
%     curlz(iz,:)=(f(iz+1,:)-f(iz-1,:))/2/dz;    %\partial/\partial z
% end %for
% curlz(lz,:)=(f(lz,:)-f(lz-1,:))/dz;
% 
% %add quiver for curl on top of color plot
% figure(3);
% contour3(x,y,f);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title('f(x,y,z) and grad(f)');
% colorbar;
% quiver3(X,Y,Z,curlx,curly,curlz,'Color','white','LineWidth',2);
% set(gca,'FontSize',24);
% hold off;

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


%z-derivative part of the divergence 
% divz=zeros(size(z));
% divz(1,:)=(g(2,:)-g(1,:))/dz;
% for iz=2:lz-1
%     divz(iz,:)=(g(iz+1,:)-g(iz-1,:))/2/dz;
% end %for
% divz(ly,:)=(g(lz,:)-g(lz-1,:))/dz;

div=divx+divy;    %this is really laplacian b/c input is gradient

figure(4);
surface(x,y,div);
set(gca,'FontSize',24);
xlabel('x');
ylabel('y');
zlabel('z');
title('laplacian(f)');
colorbar;