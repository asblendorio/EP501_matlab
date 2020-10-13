%% Multidimensional function from Problem #3
fm=@objfun2Df_AS;
gm=@objfun2Dg_AS;
gradfm=@grad_objfun2Df_AS;
gradgm=@grad_objfun2Dg_AS;

%this is for plotting
x=linspace(-1.5,1.5,20);
y=linspace(-1.5,1.5,20);
[X,Y]=meshgrid(x,y);
F=fm(X,Y);
G=gm(X,Y);

%% Newton's method for multi-variable nonlinear equations
%x0=i;
%y0=0.258*i;
x0=0.1;
y0=0.1;
[xm,ym,it2D,success2D]=newton2D_exactAS(fm,gradfm,gm,gradgm,x0,y0,100,1e-6,true);

figure;
surf(X,Y,F);
hold on;
surf(X,Y,G);
plot3(xm,ym,0,'wo','MarkerSize',32,'LineWidth',8);
