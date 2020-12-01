%% Line Integration 
% Compute and plot the parametric path 
lx=100;
ly=100;
%constants
I=10;
a=0.005;
m0=(4*pi)*10.^(-7);

x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);
[X,Y]=meshgrid(x,y);

