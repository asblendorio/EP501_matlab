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
 
%% Part A -- Handwritten 
% (a) Combine the four update steps for RK4 into a single formula an use it to derive a condition
% describing the stability of RK4. Your condition should be of the form 
% where G is a polynomial in the quantity ∆t.
disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%Part 1A Solution:%%%%%%%');

disp('%%%%%%%%End Part 1A Solution:%%%%%%%');
%% Part B
% (b) Plot your gain factor G vs. α∆t and mark or otherwise identify regions 
% of stability and instability in your plot.


disp('%%%%%%%%End Part 1B Solution:%%%%%%%');
%% Part C 
% (c) Note that the condition in Equation 7 can be expressed as two conditions (on account of the absolute value):
% −1 < G(α,∆t) < 1 (8) The marginal stability limits (the points where we transition from stable to unstable), then, are
% given by: and G(α, ∆t) − 1 = 0 (9) G(α, ∆t) + 1 = 0 (10)
% which forms a pair of root finding problems. Plot each of these marginal 
% stability conditions and evaluate all real-valued roots using an appropriate numerical method.



%% Part D
% (d) What is the largest time step for which RK4 will give stable results
% for the test ODE of Equation 6? How does that compare to RK2 stability. 



%% Part E 
% (e) Numerically solve the given ODE with RK4 using time steps slightly above
% and below your derived stability criteria and show plots for each simulation
% that demonstrate that it behaves as your analysis predicts for these two choices of time step.




disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER END%%%%%%%%%%%%%%%%%%');
%% Problem 2
disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
% Taylor series expansions were our preferred method for deriving finite difference
% approximations to derivatives. In class we derived these by repeatedly 
% generating expansions until we were able to solve for the desired derivative
% at the desired level of accuracy. There is, however, a more “organized” 
% process for deriving derivatives of arbitrary order-this problem explores one such process.
%% Part A -- Handwritten
% First and second derivative approximations were generated in class via Taylor
% expansions (see notes). Use Taylor series to generate a third derivative 
% (d^3f/dx^3) finite difference formula ap- proximation for a function f(x) 
% specified on some grid: f(xi) = {fi} having step size ∆x. There are 
% several such approximations - generate specifically one involving the 
% function values fi+2, fi+1, fi, fi−1.
%% Part B 
% Write a script to test your third derivative formula on the function:
% f(x) = cosx and plot your result alongside the analytical third derivative
% of this function (which you compute by hand). Ignore boundary adjacent 
% points (i.e. only compute the third derivative on “interior” grid points).
% Use a 100 point grid covering the region 0 ≤ x ≤ 2π.
%% Part C 
% Generally an Nth derivative requires N Taylor series expansions in order 
% to solve for the desired derivatives in terms of the gridded function values {fi}.
% For example, the fourth derivative d^4f/dx^4 can be derived from the four 
% Taylor series for fi+2,fi+1,fi−1,fi−2. Generate Taylor series for these 
% quantities in terms of the derivatives at the ith grid point (e.g. f′(xi) and so on).
%% Part D 
% Use your Taylor series equations to formulate a matrix system of equations
% for the unknowns fi′, f′′, f′′′, and f(4) (viz. the derivatives up to 
% fourth order at the ith grid point), express your system in the form:
% Where the Mjk entries are obtained from the Taylor series. In compact 
% matrix form this can be represented as: where: ∆f = M f′


%% Part E 
% Write a MATLAB or Python script to numerically invert this system using 
% either Gauss-Jordan elimination or LU factorization to compute the M−1 
% so that you have a system describing the derivatives as a function of the
% function data (differences) at various grid points:
% f′= M−1 ∆f Once this is solved for f′, one may solve for the desired
% derivative by hand as needed by dividing through by ∆x^4 and combining fi terms.
% Derive formula for the fourth derivative using your −1.

%% Part F
disp('%%%%%%%%Part 2F Solution:%%%%%%%');
% Write a function to compute the fourth derivative using the formula derived
% in part e and use it to differentiate the test function (Equation 11), 
% over interior grid points of your domain. Plot the result alongside the 
% analytical fourth derivative that you compute by hand.

disp('%%%%%%%%End Part 2F Solution:%%%%%%%');

%% Part G 
% Extend your code from part f so that it can (in principle) be used to 
% find formulas for derivatives of arbitrary order. To improve accuracy 
% use a centered stencil, e.g. for an nth order derivative use function 
% values 􏰉fi−⌈n/2⌉ ...fi ...fi+⌈n/2⌉􏰊 from grid point indices i−⌈n/2⌉...i...i+⌈n/2⌉. 
% Use your code to derive and develope formulas for the fifth and sixth derivatives
% and write these in your solution (you do not need to implement these 
% derivatives, just use your program to derive their formulae).


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
%% Part A -- Handwritten 
disp('%%%%%%%%Part 3A Solution:%%%%%%%');
% (a) The simplest BDF is just a backward Euler method where now the 
% right-hand side is evaluated at the n + 1 time level (contrast with Equation 17). 
% Such approaches are unconditionally stable for the linear test problem we are using.
% Because the right-hand side is evaluated at the n+1 time you will need to 
% solve a tridiagonal system of equations at each time step. 
% Develop, starting from equation 19 the system of equations corresponding
% to the backward Euler (in time) method.

disp('%%%%%%%%End Part 3A Solution:%%%%%%%');
%% Part B
disp('%%%%%%%%Part 3B Solution:%%%%%%%');

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
legend('exact','bwd Eul.');
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