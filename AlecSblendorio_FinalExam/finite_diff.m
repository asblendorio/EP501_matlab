%% Problem 2b -- numerical differentiation
% Analytical Plotter 
lx=100;
x=linspace(0,2*pi,lx)';
y=cos(x);
dy=-sin(x);
ddy=-cos(x);
dddy=sin(x);

figure;
plot(x,y)
hold on;
plot(x,dy)
hold on;
plot(x,ddy)
hold on;
plot(x,dddy);
legend('original function','1st derivative','2nd derivative','3rd derivative');
xlabel('x');
ylabel('y(x) or y''(x)');
title('Analytically Solved Functions');

%% Comparison of basic numerical derivative
%second order, centered
dy_dx=zeros(lx,1);
dx=x(2)-x(1);

%forward difference at the beginning
dy_dx(1)=(y(2)-y(1))/dx;

%centered difference on the interior
for ix=2:lx-1
    dy_dx(ix)=(y(ix+1)-y(ix-1))/2/dx;
end %for

%backward difference at the end
dy_dx(lx)=(y(lx)-y(lx-1))/dx;

plot(x,dy_dx,'k--')

%first order derivative approximation (backward)
%interior
dy_dxbwd1=zeros(lx,1);
dy_dxbwd2=zeros(lx,1);
dy_dxbwd3=zeros(lx,1);

for ix=2:lx
    for ij=2:lx
        for ik=2:lx
        dy_dxbwd1(ix)=(y(ix)-y(ix-1))/dx;
        dy_dxbwd2(ij)=(y(ij)-2.*y(ij-1)+y(ij-2))/2*dx;
        dy_dxbwd3(ik)=(y(ik-2)-6.*y(ik-1)+3.*y(ik)+2.*y(ik+1))/6*dx;
    
        end %for
    end %for
end %for
dy_dxbwd1(1)=dy_dxbwd1(2);
dy_dxbwd2(1)=dy_dxbwd2(2);
dy_dxbwd3(1)=dy_dxbwd3(2);

% %second order derivative approximation (backward)
% %interior
% for ix=2:lx
%     dy_dxbwd(ix)=(y(ix)-2.*y(ix-1)+y(ix-2)/2*dx;
% end %for
% dy_dxbwd(1)=dy_dxbwd(2);
% 
% %third order derivative approximation (backward)
% for ix=2:lx
%     dy_dxbwd(ix)=(y(ix-2)-6.*y(ix-1)+3.*y(ix)+2.*y(ix+1))/6*dx;
% end %for
% dy_dxbwd(1)=dy_dxbwd(2);

% plot(x,dy_dxbwd,'m--')
% legend('original function','analytical','centered','backward')
% xlabel('x');
% ylabel('y(x) or y''(x)');
% title('Comparison of finite difference derivative approximations');

plot(x,dy_dxbwd1,'m--')
hold on;
plot(x,dy_dxbwd2,'m--')
hold on;
plot(x,dy_dxbwd3,'m--')

legend('original function','analytical','centered','backward')
xlabel('x');
ylabel('y(x) or y''(x)');
title('Comparison of finite difference derivative approximations');


