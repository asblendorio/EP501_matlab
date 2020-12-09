%% Author: Alec Sblendorio 
%% Due Date: December 11, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% FINAL EXAM 
%% Import Data 

%% Problem 1 
% Problem #1: In class and in the homework we discussed the fourth order Runge-Kutta method:
% with appropriate definitions for the ∆yi as given in the book and
% repeated here for convenience
% This problem explores the stability of RK4 via several different
% approaches using the test problem.where α is a positive and real-valued constant.
 
%% Part A 
% (a) Combine the four update steps for RK4 into a single formula an use it to derive a condition
% describing the stability of RK4. Your condition should be of the form 
% where G is a polynomial in the quantity ∆t.
disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%Part 1A HANDWRITTEN Solution:%%%%%%%');

disp('%%%%%%%%End Part 1A HANDWRITTEN Solution:%%%%%%%');
%% Part B
% (b) Plot your gain factor G vs. α∆t and mark or otherwise identify regions 
% of stability and instability in your plot.
disp('%%%%%%%%Part 1B Solution:%%%%%%%');
%% Gridding in time
N=50;
tmin=0;
tmax=10;
t=linspace(tmin,tmax,N);
dt=t(2)-t(1);

%% Test problem, true solution
y0=1;
alpha=2;
ybar=y0*exp(-alpha*t);
%% RK4 Stability Considertions, FDE Analysis
adt=linspace(0.001,3.5,50);
ladt=numel(adt);
F=zeros(ladt,1);

for igain=1:ladt
    F(igain)=(1-adt(igain)+1/2*adt(igain).^2-1/6*adt(igain).^3+1/24*adt(igain).^4);
end %for

figure(2);
plot(adt,F,'*');
set(gca,'FontSize',20);
xlabel('\alpha \Delta t');
ylabel('gain factor');
title('Gain Factor vs. \alpha \Delta t');

disp('%%%%%%%%End Part 1B Solution:%%%%%%%');
%% Part C 
% (c) Note that the condition in Equation 7 can be expressed as two conditions
% (on account of the absolute value): −1 < G(α,∆t) < 1 
% The marginal stability limits (the points where we transition from stable
% to unstable), then, are given by:
% G(α, ∆t) − 1 = 0 &  G(α, ∆t) + 1 = 0 
% which forms a pair of root finding problems. Plot each of these marginal 
% stability conditions and evaluate all real-valued roots using an 
% appropriate numerical method.
disp('%%%%%%%%Part 1C Solution:%%%%%%%');
F1 = F+1;
F2 = F-1;

Fup = F1';
Flow = F2';

figure(3);
plot(adt,F-1,'*');
hold on;
plot(adt,F+1,'^');
patch([adt fliplr(adt)],[Fup fliplr(Flow)],'b--');
set(gca,'FontSize',20);
xlabel('\alpha \Delta t');
ylabel('gain factor');
title('Marginal Stability Limits and Stable Region');
legend('G-1','G+1');

%%Using the Newton Exact Method, we will evaluate all real-valued roots
%% Params for Newton iteration
maxit=100;       %maximum number of iterations allowed
minx=0;
maxx=2*pi;
tol=1e-6;        %how close to zero we need to get to cease iterations

%% Newton-Rhapson root-finding method for polynomials 
%%this method will find all of the REAL-valued roots of a polynomial
% Real roots for G - 1 = 0;
f=@objfun1;      %set the function for which we are finding roots, change to illustrate different problems
fprime=@objfun1_deriv;
verbose=true;
j=0;
finalarray=[];

for i = 0:1.0:10
    [xNewton1,itNew1,flag1]=newton_exact(f,fprime,i,maxit,tol,verbose); 
    j=j+1; 
    finalarray1(j)=xNewton1;
end

result1=finalarray1(1,2);
fprintf('The Roots of the first polynomial are: %d',result1);

disp('%%%%%%%%End Part 1C Solution:%%%%%%%');
%% Part D
% (d) What is the largest time step for which RK4 will give stable results
% for the test ODE of Equation 6? How does that compare to RK2 stability. 
disp('%%%%%%%%Part 1D Solution:%%%%%%%');
%% Gridding in time
N=1000;
tmin=0;
tmax=10;
t=linspace(tmin,tmax,N);
dt=t(2)-t(1);


%% Test problem, true solution
y0=1;
alpha=2;
ybar=y0*exp(-alpha*t);


%% Second order method; RK2
yRK2=zeros(1,N);
yRK2(1)=y0;
for n=2:N
    yhalf=yRK2(n-1)+dt/2*(-alpha*yRK2(n-1));
    yRK2(n)=yRK2(n-1)+dt*(-alpha*yhalf);
end %for


%% RK4 example; comparison against first and second order methods
yRK4=zeros(1,N);
yRK4(1)=y0;
for n=2:N
    dy1=dt*fRK(t(n-1),yRK4(n-1),alpha);
    dy2=dt*fRK(t(n-1)+dt/2,yRK4(n-1)+dy1/2,alpha);
    dy3=dt*fRK(t(n-1)+dt/2,yRK4(n-1)+dy2/2,alpha);
    dy4=dt*fRK(t(n-1)+dt,yRK4(n-1)+dy3,alpha);
    
    yRK4(n)=yRK4(n-1)+1/6*(dy1+2*dy2+2*dy3+dy4);
end %for


%% Plots of RK solutions against true solution
figure(4);
clf;
plot(t,ybar,'o-');
xlabel('t');
ylabel('y(t)');
set(gca,'FontSize',20);
figure(5);
hold on;
plot(t,yRK2,'--')
figure(1);
plot(t,yRK4,'^-')
legend('exact','RK2','RK4')

%% RK2 stability considerations, FDE analysis
adt=linspace(0.01,3,200);
ladt=numel(adt);
G=zeros(ladt,1);
for igain=1:ladt
    G(igain)=(1-adt(igain)+1/2*adt(igain).^2);
end %for
figure(6);
plot(adt,G,'o')
set(gca,'FontSize',20);
xlabel('\alpha \Delta t');
ylabel('gain factor');
title('RK2 Stability');
hold on;
%% RK4 Stability Considertions, FDE Analysis
adt2=linspace(0.01,3,200);
ladt=numel(adt2);
F=zeros(ladt,1);

for igain=1:ladt
    F(igain)=(1-adt2(igain)+1/2*adt2(igain).^2-1/6*adt2(igain).^3+1/24*adt2(igain).^4);
end %for

plot(adt2,F,'*')
set(gca,'FontSize',20);
xlabel('\alpha \Delta t');
ylabel('gain factor');
title('RK2 vs. RK4 Stability');
legend('RK2','RK4');

disp('%%%%%%%%End Part 1D Solution:%%%%%%%');
%% Part E 
% (e) Numerically solve the given ODE with RK4 using time steps slightly above
% and below your derived stability criteria and show plots for each simulation
% that demonstrate that it behaves as your analysis predicts for these two choices of time step.
disp('%%%%%%%%Part 1E Solution:%%%%%%%');


disp('%%%%%%%%End Part 1E Solution:%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER END%%%%%%%%%%%%%%%%%%');
%% Problem 2
disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
% Taylor series expansions were our preferred method for deriving finite difference
% approximations to derivatives. In class we derived these by repeatedly 
% generating expansions until we were able to solve for the desired derivative
% at the desired level of accuracy. There is, however, a more “organized” 
% process for deriving derivatives of arbitrary order - this problem explores one such process.
%% Part A 
% First and second derivative approximations were generated in class via Taylor
% expansions (see notes). Use Taylor series to generate a third derivative 
% (d3f/dx3) finite difference formula ap- proximation for a function f(x) 
% specified on some grid: f(xi) = {fi} having step size ∆x. There are 
% several such approximations - generate specifically one involving the 
% function values fi+2, fi+1, fi, fi−1.
disp('%%%%%%%%Part 2A HANDWRITTEN Solution:%%%%%%%');

disp('%%%%%%%%End Part 2A HANDWRITTEN Solution:%%%%%%%');

%% Part B 
% Write a script to test your third derivative formula on the function:
% f(x) = cosx and plot your result alongside the analytical third derivative
% of this function (which you compute by hand). Ignore boundary adjacent 
% points (i.e. only compute the third derivative on “interior” grid points).
% Use a 100 point grid covering the region 0 ≤ x ≤ 2π.
disp('%%%%%%%%Part 2B Solution:%%%%%%%');
%% Analytical Solution 
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

%% Numerical derivative
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







disp('%%%%%%%%End Part 2B Solution:%%%%%%%');
%% Part C 
% Generally an Nth derivative requires N Taylor series expansions in order 
% to solve for the desired derivatives in terms of the gridded function values {fi}.
% For example, the fourth derivative d4f/dx4 can be derived from the four 
% Taylor series for fi+2,fi+1,fi−1,fi−2. Generate Taylor series for these 
% quantities in terms of the derivatives at the ith grid point (e.g. f′(xi)
% and so on).
disp('%%%%%%%%Part 2C Solution:%%%%%%%');

disp('%%%%%%%%End Part 2C Solution:%%%%%%%');
%% Part D 
% Use your Taylor series equations to formulate a matrix system of equations
% for the unknowns fi′, f′′, f′′′, and f(4) (viz. the derivatives up to 
% fourth order at the ith grid point), express your system in the form:
% Where the Mjk entries are obtained from the Taylor series. In compact 
% matrix form this can be represented as: where: ∆f = M f′
disp('%%%%%%%%Part 2D Solution:%%%%%%%');

disp('%%%%%%%%End Part 2D Solution:%%%%%%%');

%% Part E 
% Write a MATLAB or Python script to numerically invert this system using 
% either Gauss-Jordan elimination or LU factorization to compute the M−1 
% so that you have a system describing the derivatives as a function of the
% function data (differences) at various grid points:
% f′ =M−1 ∆f (15) Once this is solved for f′, one may solve for the desired
% derivative by hand as needed by dividing through by ∆x4 and combining fi terms.
% Derive formula for the fourth derivative using your −1.
disp('%%%%%%%%Part 2E Solution:%%%%%%%');

disp('%%%%%%%%End Part 2E Solution:%%%%%%%');
%% Part F
% Write a function to compute the fourth derivative using the formula derived
% in part e and use it to differentiate the test function (Equation 11), 
% over interior grid points of your domain. Plot the result alongside the 
% analytical fourth derivative that you compute by hand.
disp('%%%%%%%%Part 2F Solution:%%%%%%%');

disp('%%%%%%%%End Part 2F Solution:%%%%%%%');

%% Part G 
% Extend your code from part f so that it can (in principle) be used to 
% find formulas for derivatives of arbitrary order. To improve accuracy 
% use a centered stencil, e.g. for an nth order derivative use function 
% values 􏰉fi−⌈n/2⌉ ...fi ...fi+⌈n/2⌉􏰊 from grid point indices i−⌈n/2⌉...i...i+⌈n/2⌉. 
% Use your code to derive and develope formulas for the fifth and sixth derivatives
% and write these in your solution (you do not need to implement these 
% derivatives, just use your program to derive their formulae).
disp('%%%%%%%%Part 2G Solution:%%%%%%%');

disp('%%%%%%%%End Part 2G Solution:%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER END%%%%%%%%%%%%%%%%%%');
%% Problem 3 
% Parabolic partial differential equations of the form are often difficult 
% to deal with due to the highly restrictive stability condition applying 
% to forward difference (in time) methods like a forward Euler approach 
% (note the right-hand side is evaluated at the n time level)
% In practice methods based on forward differencing in time are often 
% unstable for reasonable time steps. Irrespective of the time differencing
% used, the spatial derivative can be any reasonable finite difference
% approximation; for this problem we use a second order accurate centered 
% difference (here specified at an arbitrary mth time level):
% Because of overly restrictive stability criteria for forward in time differences, 
% backward time difference formulas (BDFs) are more commonly used for parabolic equations.
% This problem explores several commonly-used BDFs that will be applied to
% solve Equation 16.

disp('%%%%%%%%%%%%%%%%%%PROBLEM #3 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
%% Part A 
% (a) The simplest BDF is just a backward Euler method where now the 
% right-hand side is evaluated at the n + 1 time level (contrast with Equation 17). 
% Such approaches are unconditionally stable for the linear test problem we are using.
% Because the right-hand side is evaluated at the n+1 time you will need to 
% solve a tridiagonal system of equations at each time step. 
% Develop, starting from equation 19 the system of equations corresponding
% to the backward Euler (in time) method.
disp('%%%%%%%%Part 3A Solution:%%%%%%%');

disp('%%%%%%%%End Part 3A Solution:%%%%%%%');
%% Part B
% Write a MATLAB or Python script that solves Equation 16 using a backward
% Euler in time method with a centered in space derivatives (Equation 18). 
% You may use built-in MATLAB or Python utilities for the matrix solution 
% or you can leverage your tridiagonal solver from the midterm or the 
% Gaussian elimination tools from the repositories. Note that both repositories
% have an example of the Crank-Nicolson method (in MATLAB or Python) 
% which is quite similar to backward Euler and may serve as a useful template. 
% Plot your results in a manner similar to what is done for the example script provided.
% Use the parameters: λ=2; ∆x = 1/64; ∆t = 5*(∆x^2/2λ); 0≤x≤1 ; 0≤t≤1024
% 1/λ(2π(2∆x)^2)

disp('%%%%%%%%Part 3B Solution:%%%%%%%');

%% Gridding in time
N=25;
tmin=0;
tmax=6;
t=linspace(tmin,tmax,N);
dt=t(2)-t(1);


%% Test problem analytical solution
y0=1;
alpha=2;
ybar=y0*exp(-alpha*t);

%% Backward Euler solution
ybwd=zeros(1,N);
ybwd(1)=y0;
for n=2:N
    ybwd(n)=ybwd(n-1)/(1+alpha*dt);
end %for


%% Plot results for all solutions
figure(1);
plot(t,ybar,'o-');
xlabel('t');
ylabel('y(t)');
set(gca,'FontSize',20);
hold on;
plot(t,ybwd,'-.');
legend('exact','bwd Eul.')

disp('%%%%%%%%End Part 3B Solution:%%%%%%%');
%% Part C
% The backward Euler method is only first order accurate in time. 
% A second order in time BDF approach can be developed using a second order
% accurate backward difference for the time derivative.
% Develop a numerical solution to Equation 16 based on this approach. 
% Note that you will need to store data for fi at two previous time levels 
% (n and n − 1) as opposed to the backward Euler approach which only requires
% prior data from the n time level. Develop and write down your system of 
% equations that would need to be solved at each time step.
disp('%%%%%%%%Part 3C Solution:%%%%%%%');

disp('%%%%%%%%End Part 3C Solution:%%%%%%%');

%%  Part D  
% Make a version of your parabolic solver from part b that implements the 
% second order BDF derived in part c and run it on the same test problem. 
% Plot your results and compare them against the backward Euler approach 
% and analytical solution.
disp('%%%%%%%%Part 3D Solution:%%%%%%%');

disp('%%%%%%%%End Part 3D Solution:%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%PROBLEM #3 ANSWER END%%%%%%%%%%%%%%%%%%');