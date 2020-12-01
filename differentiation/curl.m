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
hold off;

figure(2);
contourf(x,y,f);
xlabel('x');
ylabel('y');
xlim([-0.02 0.02]);
ylim([-0.02 0.02]);
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

%forward difference at the beginning
dy_dx(1)=(y(2)-y(1))/dx;

%centered difference on the interior
for ix=2:lx-1
    dy_dx(ix)=(y(ix+1)-y(ix-1))/2/dx;
end %for

%backward difference at the end
dy_dx(lx)=(y(lx)-y(lx-1))/dx;

plot(x,dy_dx,'k--');

curl_B=curlx-curly;
figure(3);
imagesc(x,y,curl_B);
xlabel('x');
ylabel('y');
title('Curl of B');
hold on;
colorbar;
%quiver(X,Y,dx_dy,dy_dx,'Color','white','LineWidth',2);
set(gca,'FontSize',24);
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

