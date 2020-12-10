%% Problem 2b 
% Analytical Plotter 
lx=100;
tmin=0;
tmax=2*pi;
x=linspace(tmin,tmax,lx);
dx=x(2)-x(1);
y=cos(x);
dy=-sin(x);
ddy=-cos(x);
dddy=sin(x);

figure(1);
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

%% numerical derivative
%first,second,third order derivative approximation (backward)
%interior
dx1=zeros(lx,1);
dx2=zeros(lx,1);
dx3=zeros(lx,1);
dx4=zeros(lx,1);

for ix=3:lx-2
        dx1(ix)=(y(ix)-y(ix-1))/dx;   
end %for
dx1(1)=dx1(2);

for ix=3:lx-2
    dx2(ix)=(y(ix)-2.*y(ix-1)+y(ix-2))/2*dx;
end %for
dx2(1)=dx2(2);

for ix=3:lx-2
    dx3(ix)=(y(ix+2)-6.*y(ix-1)+3.*y(ix)+2.*y(ix+1))/6*dx;
end %for    
dx3(1)=dx3(2);

for ix=3:lx-2
    dx4(ix)=(y(ix-2)-8.*y(ix-1)+8.*y(ix+1)-y(ix+2))/12*dx;
end %for    
dx4(1)=dx4(2);

figure(2);
plot(x,dx1,'m--')
hold on;
plot(x,dx2,'b--')
hold on;
plot(x,dx3,'k--')
hold on;
plot(x,dx4,'r--')

legend('1st derivative','2nd derivative','3rd derivative','4th derivative');
xlabel('x');
ylabel('y(x) or y''(x)');
title('Numerically Solved Functions');

