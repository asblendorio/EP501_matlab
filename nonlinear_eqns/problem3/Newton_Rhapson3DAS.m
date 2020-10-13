%% Multidimensional function from Problem #3
%%Part 3b of problem, compute the roots of 3 equations (3 dimensions)
fm=@objfun3Df_AS;
gm=@objfun3Dg_AS;
km=@objfun3Dk_AS;
gradfm=@grad_objfun3Df_AS;
gradgm=@grad_objfun3Dg_AS;
gradkm=@grad_objfun3Dk_AS;
%this is for plotting
x=linspace(-1.5,1.5,20);
y=linspace(-1.5,1.5,20);
z=linspace(-1.5,1.5,20);
[X,Y,Z]=meshgrid(x,y,z);
F=fm(X,Y,Z);
G=gm(X,Y,Z);
K=km(X,Y,Z);

%% Newton's method for multi-variable nonlinear equations
%x0=i;
%y0=0.258*i;
%z0=0.58*i;
x0=0.1;
y0=0.1;
z0=0.1;
[xm,ym,zm,it3D,success3D]=newton3D_exactAS(fm,gradfm,gm,gradgm,km,gradkm,x0,y0,z0,100,1e-6,true);

figure;
surf(X,Y,F);
hold on;
surf(X,Y,G);
hold on;
surf(X,Y,K);
plot3(xm,ym,0,'wo','MarkerSize',32,'LineWidth',8);
