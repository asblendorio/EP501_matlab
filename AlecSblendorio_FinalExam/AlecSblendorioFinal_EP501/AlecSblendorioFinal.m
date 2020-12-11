%% Author: Alec Sblendorio 
%% Due Date: December 11, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% FINAL EXAM 
%% Import Data 
load('matrix_M.mat');
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
f = F-1;
figure(2);
plot(adt,f,'*');
hold on;
x1=0:4;
y1=0;
plot(x1,y1*ones(size(x1)),'LineWidth',3);
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
F1 = f+1;
F2 = f-1;

Fup = F1';
Flow = F2';

figure(3);
plot(adt,F2,'*');
hold on;
plot(adt,F1,'^');
hold on;
%Plotting horizontal line to indicate limit of stability
x2=0:4;
y2=0;
plot(x2,y2*ones(size(x2)),'LineWidth',3);

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

for i = 0:1.0:10
    [xNewton,itNew,flag]=newton_exact(f,fprime,i,maxit,tol,verbose); 
    j=j+1; 
    finalarray1(j)=xNewton;
end

result1=0;
result2=finalarray1(1,2);
disp('Please Ignore the Warning from Newton Exact, the program seems to be working well.');
fprintf('The Real Roots of the polynomial are: %d and %f\n',result1,result2);

disp('%%%%%%%%End Part 1C Solution:%%%%%%%');
%% Part D
% (d) What is the largest time step for which RK4 will give stable results
% for the test ODE of Equation 6? How does that compare to RK2 stability. 
disp('%%%%%%%%Part 1D Solution:%%%%%%%');
%% Gridding in time
N=250;
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


%% RK2 stability considerations, FDE analysis
adt=linspace(0.01,3,200);
ladt=numel(adt);
G=zeros(ladt,1);
for igain=1:ladt
    G(igain)=(1-adt(igain)+1/2*adt(igain).^2);
end %for
figure(4);
plot(adt,G-1,'o')
set(gca,'FontSize',20);
hold on;
%% RK4 Stability Considertions, FDE Analysis
adt2=linspace(0.01,3,200);
ladt=numel(adt2);
F=zeros(ladt,1);

for igain=1:ladt
    F(igain)=(1-adt2(igain)+1/2*adt2(igain).^2-1/6*adt2(igain).^3+1/24*adt2(igain).^4);
end %for

plot(adt2,F-1,'*')
hold on;
x3=0:4;
y3=0;
plot(x3,y3*ones(size(x3)),'LineWidth',3);
set(gca,'FontSize',20);
xlabel('\alpha \Delta t');
ylabel('gain factor');
title('RK2 vs. RK4 Stability Curves');
legend('RK2','RK4');
hold off;
disp('As seen in Figure 4, the largest time step for which RK4 will give stable results for the test ODE is 2.785');
disp('RK4 has a higher time step than RK2 showing the 4th order method will contain a tighter tolerance when evaluating a function. If we evaluate the same function with the 2nd order method, it will not produce similar tolerance.');

%% Gridding in time for RK4
N_4=8;
tmin=0;
tmax_4=10;
t4=linspace(tmin,tmax_4,N_4);
dt_4=t4(2)-t4(1);

%% Gridding in time for RK2
N_2=8;
tmin=0;
tmax_2=10;
t2=linspace(tmin,tmax_2,N_2);
dt_2=t2(2)-t2(1);
%% Test problem, true solution Rk4
y0=1;
alpha=1;
ybar4=y0*exp(-alpha*t4);
%% Test problem, true solution Rk2
ybar2=y0*exp(-alpha*t2);

%% Second order method; RK2
yRK2=zeros(1,N_2);
yRK2(1)=y0;
for n=2:N_2
    yhalf=yRK2(n-1)+dt_2/2*(-alpha*yRK2(n-1));
    yRK2(n)=yRK2(n-1)+dt_2*(-alpha*yhalf);
end %for
%% RK4 example; comparison against first and second order methods
yRK4=zeros(1,N_4);
yRK4(1)=y0;
for n=2:N_4
    dy1=dt_4*fRK(t(n-1),yRK4(n-1),alpha);
    dy2=dt_4*fRK(t(n-1)+dt_4/2,yRK4(n-1)+dy1/2,alpha);
    dy3=dt_4*fRK(t(n-1)+dt_4/2,yRK4(n-1)+dy2/2,alpha);
    dy4=dt_4*fRK(t(n-1)+dt_4,yRK4(n-1)+dy3,alpha);
    
    yRK4(n)=yRK4(n-1)+1/6*(dy1+2*dy2+2*dy3+dy4);
end %for

%% Plots of RK solutions against true solution
figure(6);
clf;
plot(t4,ybar4,'o-');
xlabel('t');
ylabel('y(t)');
set(gca,'FontSize',20);

figure(6);
hold on;
plot(t4,yRK4,'^-')
title('Time Step Analysis for RK4');
legend('exact','RK4 N=8')

figure(7);
plot(t2,ybar2,'o-');
xlabel('t');
ylabel('y(t)');
set(gca,'FontSize',20);
figure(7);
hold on;
plot(t2,yRK2,'--')
title('Time Step Analysis for RK2');
legend('exact','RK2 N=8')

fprintf('The largest time step for RK4 to maintain stability is: %d.\n',N_4);
fprintf('the largest time step for RK2 to maintain stability is: %d\n',N_2);
%% Comments on 1D 
% In Figures 6 and 7, I present both methods with the same time step. It is
% clear that at that specific time step, the RK 4 method has greater
% stability than RK2. 
disp('%%%%%%%%End Part 1D Solution:%%%%%%%');
%% Part E 
% (e) Numerically solve the given ODE with RK4 using time steps slightly above
% and below your derived stability criteria and show plots for each simulation
% that demonstrate that it behaves as your analysis predicts for these two choices of time step.
disp('%%%%%%%%Part 1E Solution:%%%%%%%');
%% Gridding in time for RK4 Above derived stability Criteria
N_4above=80;
tmin=0;
tmax_4above=10;
t4above=linspace(tmin,tmax_4above,N_4above);
dt_4above=t4above(2)-t4above(1);
%% Gridding in time for RK4 Below derived stability Criteria
N_4below=4;
tmin=0;
tmax_4below=10;
t4below=linspace(tmin,tmax_4below,N_4below);
dt_4below=t4below(2)-t4below(1);
%% Test problem, true solution Rk4
y0=1;
alpha=1;
ybar4_above=y0*exp(-alpha*t4above);
ybar4_below=y0*exp(-alpha*t4below);

%% RK4 Above  example; comparison against first and second order methods
yRK4_above=zeros(1,N_4above);
yRK4_above(1)=y0;
for n=2:N_4above
    dy1=dt_4above*fRK(t(n-1),yRK4_above(n-1),alpha);
    dy2=dt_4above*fRK(t(n-1)+dt_4/2,yRK4_above(n-1)+dy1/2,alpha);
    dy3=dt_4above*fRK(t(n-1)+dt_4/2,yRK4_above(n-1)+dy2/2,alpha);
    dy4=dt_4above*fRK(t(n-1)+dt_4,yRK4_above(n-1)+dy3,alpha);
    
    yRK4_above(n)=yRK4_above(n-1)+1/6*(dy1+2*dy2+2*dy3+dy4);
end %for

%% RK4 Below example; comparison against first and second order methods
yRK4_below=zeros(1,N_4below);
yRK4_below(1)=y0;
for n=2:N_4below
    dy1=dt_4below*fRK(t(n-1),yRK4_below(n-1),alpha);
    dy2=dt_4below*fRK(t(n-1)+dt_4below/2,yRK4_below(n-1)+dy1/2,alpha);
    dy3=dt_4below*fRK(t(n-1)+dt_4below/2,yRK4_below(n-1)+dy2/2,alpha);
    dy4=dt_4below*fRK(t(n-1)+dt_4below,yRK4_below(n-1)+dy3,alpha);
    
    yRK4_below(n)=yRK4_below(n-1)+1/6*(dy1+2*dy2+2*dy3+dy4);
end %for
%% Plots of RK solutions against true solution
figure(8);
clf;
plot(t4above,ybar4_above,'o-');
hold on;
xlabel('t');
ylabel('y(t)');
set(gca,'FontSize',20);
plot(t4above,yRK4_above,'^-')
title('Altered Time Step Analysis for High N Value -- RK4');
legend('exact','RK4 N=80')

figure(9);
clf;
plot(t4below,ybar4_below,'o-');
xlabel('t');
ylabel('y(t)');
set(gca,'FontSize',20);
hold on;
plot(t4below,yRK4_below,'^-')
title('Altered Time Step Analysis for Low N Value -- RK4');
legend('exact','RK4 N=4')

%% Comments on 1E 
% For a low N value, the solution will destabilize. For a large N it will 
% converge to the true solution. 
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
% f(x) = cos(x) and plot your result alongside the analytical third derivative
% of this function (which you compute by hand). Ignore boundary adjacent 
% points (i.e. only compute the third derivative on “interior” grid points).
% Use a 100 point grid covering the region 0 ≤ x ≤ 2π.
disp('%%%%%%%%Part 2B Solution:%%%%%%%');
%% Analytical Solution 
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

figure(10);
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
lx1=100;
tmin=0;
tmax=2*pi;
x1=linspace(tmin,tmax,lx1);
dx_new=x1(2)-x1(1);
dx1=zeros(lx1,1);
dx2=zeros(lx1,1);
dx3=zeros(lx1,1);
dx4=zeros(lx1,1);

for ix=3:lx-2
        dx1(ix)=(y(ix)-y(ix-1))/(dx_new);
        dx2(ix)=(y(ix+1)-2*y(ix)+y(ix-1))/(dx_new^2);
        dx3(ix)=(2*y(ix-1)+3*y(ix)-6*y(ix+1)+y(ix+2))/(6*dx_new);
        dx4(ix)=(y(ix-2)-8*y(ix-1)+8*y(ix+1)-y(ix+2))/(12*dx_new);
end %for


figure(11);
plot(x1,y,'y--')
hold on;
plot(x1,dx1,'m--')
hold on;
plot(x1,dx2,'b--')
hold on;
plot(x1,dx3,'k--')
hold on;
plot(x,dx4,'g--')

legend('Original Function','1st derivative','2nd derivative','3rd derivative','4th Derivative');
xlabel('x');
ylabel('y(x) or y''(x)');
title('Numerically Solved Functions');

%% Comments on 2B 
% There is a steep line on the left hand and right hand side of the
% derivative functions. I don't know how to get rid of it but the rest of
% the numerical functions match the analytical. 
disp('%%%%%%%%End Part 2B Solution:%%%%%%%');
%% Part C 
% Generally an Nth derivative requires N Taylor series expansions in order 
% to solve for the desired derivatives in terms of the gridded function
% values. For example, the fourth derivative d4f/dx4 can be derived from the four 
% Taylor series for fi+2,fi+1,fi−1,fi−2. Generate Taylor series for these 
% quantities in terms of the derivatives at the ith grid point.
disp('%%%%%%%%Part 2C HANDWRITTEN Solution:%%%%%%%');

disp('%%%%%%%%End Part 2C HANDWRITTEN Solution:%%%%%%%');
%% Part D 
% Use your Taylor series equations to formulate a matrix system of equations
% for the unknowns fi′, f′′, f′′′, and f(4) (viz. the derivatives up to 
% fourth order at the ith grid point), express your system in the form:
% Where the Mjk entries are obtained from the Taylor series. In compact 
% matrix form this can be represented as: where: ∆f = M f′
disp('%%%%%%%%Part 2D HANDWRITTEN Solution:%%%%%%%');

disp('%%%%%%%%End Part 2D HANDWRITTEN Solution:%%%%%%%');

%% Part E 
% Write a MATLAB or Python script to numerically invert this system using 
% either Gauss-Jordan elimination or LU factorization to compute the M−1 
% so that you have a system describing the derivatives as a function of the
% function data (differences) at various grid points:
% f′ =M−1 ∆f (15) Once this is solved for f′, one may solve for the desired
% derivative by hand as needed by dividing through by ∆x4 and combining fi terms.
% Derive formula for the fourth derivative using your inverse.
disp('%%%%%%%%Part 2E NUMERICAL and HANDWRITTEN Solution:%%%%%%%');
soln=@DoolittleLU;
soln1=soln(M);

%% Comments on Part 2E 
% I was not able to properly calculate the inverse of M however, this did
% not stop me from deriving a formula for the 4th derivative. Due to the
% fact that my numerically computed inverse is not fully correct, I did
% have to utilize the built-in Matlab function to 1) Validate my inverse
% formula and 2) Take the coefficients to solve the rest of the problem. 
% The bottom(4th) row in the inverse matrix are the coefficients to the derivative formula.  

disp('%%%%%%%%End Part 2E NUMERICAL and HANDWRITTEN Solution:%%%%%%%');
%% Part F
% Write a function to compute the fourth derivative using the formula derived
% in part e and use it to differentiate the test function (Equation 11), 
% over interior grid points of your domain. Plot the result alongside the 
% analytical fourth derivative that you compute by hand.
disp('%%%%%%%%Part 2F Solution:%%%%%%%');

lx4=100;
tmin4=0;
tmax4=2*pi;
x4=linspace(tmin4,tmax4,lx4);
y4=cos(x4);
dx_new=x4(2)-x4(1);
dx4th=zeros(lx4,1);

for ix=3:lx4-2
    dx4th(ix)=(y4(ix-2)-4.*y4(ix-1)-4.*y4(ix+1)+y4(ix+2)+6.*y4(ix))/(dx_new^4);
end %for

%for problem 2f
figure(15);
plot(x4,dx4th,'k--','LineWidth',1.5);
hold on;
plot(x4,y4,'r--');
legend('Numerical 4th Dv','Analytical 4th Dv');
xlabel('x');
ylabel('y(x) or y''(x)');
% xlim([0 6]);
% ylim([-1 1]);
set(gca,'FontSize',20);
title('Analytically vs. Numerically Solved 4th Derivative');

disp('%%%%%%%%End Part 2F Solution:%%%%%%%');

%% Part G 
% Extend your code from part f so that it can (in principle) be used to 
% find formulas for derivatives of arbitrary order. To improve accuracy 
% use a centered stencil. Use your code to derive and develope formulas for
% the fifth and sixth derivatives and write these in your solution 
% (you do not need to implement these derivatives, just use your program to
% derive their formulae).
disp('%%%%%%%%Part 2G Solution:%%%%%%%');
disp('Solution not found prior to deadline.');
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
disp('%%%%%%%%Part 3A HANDWRITTEN Solution:%%%%%%%');

disp('%%%%%%%%End Part 3A HANDWRITTEN Solution:%%%%%%%');
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

%% Define a 1D space and time grid in x for test problem
lx=64;
a=0;     %here a,b are the endpoints of the x-domain
b=1;     %use a square region for a test problem
x=linspace(a,b,lx);
dx=x(2)-x(1);        %grid spacing

%% Define parameters of the parabolic equation, time variable
lambda=2;
tau=1/(2*pi/(2*dx))^2/lambda;    %diffusion time scale for the equation, based on smallest resolvable spatial mode
%dt=tau/5;              %time step

dtmargin=(5/2).*(dx^2/lambda);
dt=50*dtmargin;
tmin=0;
tmax=1024*tau;          %go out to three times the diffusion time scale for the smallest possible mode
t=tmin:dt:tmax;
lt=numel(t);

%% FTCS implementation
f=zeros(lx,lt);
f(:,1)=sin(2*pi*x)+sin(8*pi*x);

%FTCS iterations
for n=1:lt-1
    f(1,n+1)=0;   %assume temperature goes to some small number on left boundary
    for i=2:lx-1     %interior grid points
        f(i,n+1)=f(i,n)+dt/dx^2*lambda*(f(i+1,n)-2*f(i,n)+f(i-1,n));
    end %for
    f(lx,n+1)=0;  %assume temperature goes to some small number on right boundary
end %for

figure(13);
subplot(131);
imagesc(t,x,f);
colorbar;
axis xy;
xlabel('time (s)');
ylabel('x (m)')
title('FTCS')
set(gca,'FontSize',16);

%% Trapezoidal implementation, note matrix solutions are more efficiently handled thru tri-diagonal solver; Matlab built-in will detect automatically
f2=zeros(lx,lt);
A=sparse(lx,lx);   %allocate sparse array storage (this matrix is to be tridiag)
b=zeros(lx,1);

f2(:,1)=sin(2*pi*x)+sin(8*pi*x);
for n=2:lt-1
    A(1,1)=1;
    b(1)=0;
    for ix=2:lx-1
        %i-1 coeff
        A(ix,ix-1)=-lambda/dx^2;
        
        %i coeff
        A(ix,ix)=1/dt-2*lambda/dx^2;
        
        %i+1 coeff
        A(ix,ix+1)=-lambda/dx^2;
        
        b(ix)=f2(ix,n-1)/dt+(f2(ix+1,n-1)-2*f2(ix,n-1)+f2(ix-1,n-1))/dx^2*(lambda/2);
            
    end %for
    A(lx,lx)=1;
    b(lx)=0;
    
    fnow=A\b;
    f2(:,n)=fnow;
end %for

%% Compare two solutions on plot
figure(13);
subplot(132);
imagesc(t,x,f2);
colorbar;
axis xy;
xlabel('time (s)');
ylabel('x (m)');
title('C-N solution');
set(gca,'FontSize',16);


%% Compute and plot the analytical solution (see course repository ./test_problems/ for derivation)
[T,X]=meshgrid(t,x);
tempexact=exp(-4*pi^2*lambda*T).*sin(2*pi*X)+exp(-64*pi^2*lambda*T).*sin(8*pi*X);

figure(13);
subplot(133);
imagesc(t,x,tempexact);
colorbar;
axis xy;
xlabel('time (s)');
ylabel('x (m)');
title('Exact');
set(gca,'FontSize',16);

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
disp('%%%%%%%%Part 3C NUMERICAL and HANDWRITTEN Solution:%%%%%%%');
%% Comments on 3C
% After calculating the handwritten component to solve for the previous two time levels
% (n and n-1), I ran into quite difficulty developing a numerical method to
% solve those system of equations. Below is a list of several avenues which
% I went down (theoretically) to potentially solve this problem. 
% 1) Since we get a second order in time BDF, which has 4 coefficients, the 
% backward Euler method will not work. There are two for the i index and then
% i+1 and i-1 respectively. Initially, I though that I could solve the two 
% time levels independently. However, I then realized that we
% do not know our n+1 values. Cannot solve an equation with more than two
% unknowns! 
% 2) Since our tridiagonal solver is solving a system of equations each
% iteration, and we have dependent variables, another option would be to
% concatenate two A matrices. Therefore, including all the necessary
% coefficients for our system of equations. The problem I ran into is how to
% properly edit the trapezoidal implementation. 
% 3) Another possible avenue that I tried is adjusting the variable b or
% "b2" here. In part 3b, that variable stored the n-1 time level but in
% this problem we should be able to input n+1 and then n. However, due to
% the restraints on the system of equations the dependency of
% the space and time functions does not give me an avenue to properly solve
% this problem. 

%% Define a 1D space and time grid in x for test problem
lx=64;
a=0;     %here a,b are the endpoints of the x-domain
b=1;     %use a square region for a test problem
x=linspace(a,b,lx);
dx=x(2)-x(1);        %grid spacing

%% Define parameters of the parabolic equation, time variable
lambda=2;
tau=1/(2*pi/(2*dx))^2/lambda;    %diffusion time scale for the equation, based on smallest resolvable spatial mode
%dt=tau/5;              %time step

dtmargin1=(5/2).*(dx^2/lambda);
dt1=0.5*dtmargin1;
tmin=0;
tmax=1024*tau;          %go out to three times the diffusion time scale for the smallest possible mode
t=tmin:dt:tmax;
lt=numel(t);

%% FTCS implementation
f=zeros(lx,lt);
f(:,1)=sin(2*pi*x)+sin(8*pi*x);

%FTCS iterations
for n=1:lt-1
    f(1,n+1)=0;   %assume temperature goes to some small number on left boundary
    for i=2:lx-1     %interior grid points
        f(i,n+1)=f(i,n)+dt/dx^2*lambda*(f(i+1,n)-2*f(i,n)+f(i-1,n));
    end %for
    f(lx,n+1)=0;  %assume temperature goes to some small number on right boundary
end %for

figure(14);
subplot(131);
imagesc(t,x,f);
colorbar;
axis xy;
xlabel('time (s)');
ylabel('x (m)')
title('Part 3C: FTCS')
set(gca,'FontSize',16);

%% N-1:Trapezoidal implementation, note matrix solutions are more efficiently handled thru tri-diagonal solver; Matlab built-in will detect automatically
f3=zeros(lx,lt);
A2=sparse(lx,lx);   %allocate sparse array storage (this matrix is to be tridiag)
b2=zeros(lx,1);

f3(:,1)=sin(2*pi*x)+sin(8*pi*x);
for n=2:lt-1
    A2(1,1)=1;
    b2(1)=0;
    for ix=2:lx-1
        %i-1 coeff
        A2(ix,ix-1)=lambda/dx^2;
        
        %i coeff
        A2(ix,ix)=(3/(2*dt)-2*lambda/dx^2);
        
        %i+1 coeff
        A2(ix,ix+1)=lambda/dx^2;
        
        b2(ix)=3*f3(ix,n+1)/dt+(f3(ix+1,n-1)-2*f3(ix,n-1)+f3(ix-1,n-1))/dx^2;
            
    end %for
    A2(lx,lx)=1;
    b2(lx)=0;
    
    fnow2=A2\b2;
    f3(:,n)=fnow2;
end %for

%% Compare two solutions on plot
figure(14);
subplot(132);
imagesc(t,x,f2);
colorbar;
axis xy;
xlabel('time (s)');
ylabel('x (m)');
title('Part 3C: C-N solution');
set(gca,'FontSize',16);


%% Compute and plot the analytical solution (see course repository ./test_problems/ for derivation)
[T,X]=meshgrid(t,x);
tempexact=exp(-4*pi^2*lambda*T).*sin(2*pi*X)+exp(-64*pi^2*lambda*T).*sin(8*pi*X);

figure(14);
subplot(133);
imagesc(t,x,tempexact);
colorbar;
axis xy;
xlabel('time (s)');
ylabel('x (m)');
title('Part 3C: Exact');
set(gca,'FontSize',16);

disp('%%%%%%%%End Part 3C NUMERICAL and HANDWRITTEN Solution:%%%%%%%');

%%  Part D  
% Make a version of your parabolic solver from part b that implements the 
% second order BDF derived in part c and run it on the same test problem. 
% Plot your results and compare them against the backward Euler approach 
% and analytical solution.
disp('%%%%%%%%Part 3D Solution:%%%%%%%');
disp('Solution not found prior to deadline.');
disp('%%%%%%%%End Part 3D Solution:%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%PROBLEM #3 ANSWER END%%%%%%%%%%%%%%%%%%');