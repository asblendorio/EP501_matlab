%% Problem 2b -- numerical differentiation
% Analytical Plotter 
lx=100;
x=linspace(0,2*pi,lx)';
dx=x(2)-x(1);
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
%first,second,third order derivative approximation (backward)
%interior
dx1=zeros(lx,1);
dx2=zeros(lx,1);
dx3=zeros(lx,1);

for ix=2:lx
        dx1(ix)=(y(ix)-y(ix-1))/dx;
        dx2(ix)=(y(ix)-2.*y(ix-1)+y(ix-2))/2/dx;
        dx3(ix)=(y(ix-2)-6.*y(ix-1)+3.*y(ix)+2.*y(ix+1))/6*dx;
end %for
dx1(1)=dx1(2);
dx2(1)=dx2(2);
dx3(1)=dx3(2);

figure(2);
plot(x,dx1,'m--')
hold on;
plot(x,dx2,'b--')
hold on;
plot(x,dx3,'k--')

legend('original function','analytical','centered','backward')
xlabel('x');
ylabel('y(x) or y''(x)');
title('Comparison of finite difference derivative approximations');


