%% Illustration of bilinear interpolation, single interval of interest
x=xg;
y=yg;
f=f2D;
[X,Y] = meshgrid(x,y);
[yi,y,yi1,xi,x,xi1] = bilinear(2,yg,ygi,xg,xgi);
x1=xprime;
y1=yprime;

% Manually written
fvec=f(:);
xvec=X(:);
yvec=Y(:);
M=[ones(4,1),xvec(:),yvec(:),xvec(:).*yvec(:)];
[Mmod,order]=Gauss_elim(M,fvec);
avec=backsub(Mmod(order,:));
finterpmanual=avec(1)+avec(2)*x1+avec(3)*y1+avec(4)*x1*y1;

% Matlab version
finterp=interp2(X,Y,f,x1,y1);

%Visual interpolation point:  x1,y1
figure(1);
imagesc(x,y,f);
axis xy;
xlabel('x');
ylabel('y');
c=colorbar;
ylabel(c,'f(x,y)')
hold on;
plot(xvec,yvec,'w^','MarkerSize',15,'MarkerFaceColor','white');
plot(x1,y1,'wo','MarkerSize',20,'MarkerFaceColor','white');
hold off;


%% Illustrate cubic spline approximations using Matlab functions
x=linspace(-5,5,15);
y=sin(x);
figure(2);
plot(x,y,'o','MarkerSize',20);

splinedef=spline(x,y);
x2=linspace(min(x),max(x),256);
y2=ppval(splinedef,x2);
hold on;
plot(x2,y2,'.');

y2true=sin(x2);
plot(x2,y2true,'-');

