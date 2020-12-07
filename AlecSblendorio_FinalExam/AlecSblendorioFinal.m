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
% process for deriving derivatives of arbitrary order - this problem explores one such process.
%% Part A 
% First and second derivative approximations were generated in class via Taylor
% expansions (see notes). Use Taylor series to generate a third derivative 
% (d3f/dx3) finite difference formula ap- proximation for a function f(x) 
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
% For example, the fourth derivative d4f/dx4 can be derived from the four 
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
% f′ =M−1 ∆f (15) Once this is solved for f′, one may solve for the desired
% derivative by hand as needed by dividing through by ∆x4 and combining fi terms.
% Derive formula for the fourth derivative using your −1.

%% Part F
disp('%%%%%%%%Part 2F Solution:%%%%%%%');
% Write a function to compute the fourth derivative using the formula derived
% in part e and use it to differentiate the test function (Equation 11), 
% over interior grid points of your domain. Plot the result alongside the 
% analytical fourth derivative that you compute by hand.

%% Params for Newton iteration
maxit=100;       %maximum number of iterations allowed
minx=0;
maxx=pi/4;
tol=1e-9;        %how close to zero we need to get to cease iterations

%% Objective function defs.
f=@objfun;      %set the function for which we are finding roots, change to illustrate different problems
fprime=@objfun_deriv;
x=linspace(minx,maxx,64);   %grid for basic plotting purposes
ygrid=f(x);
verbose=true;

%% Newton-Rhapson root-finding method
j=0;
rec = 0;
finalarray=[];
verbose=true;

for i = -10:0.15:10
    [xNewton,itNew,flag]=newton_exact(f,fprime,i,100,tol,verbose);    
    j=j+1; 
    finalarray(j)=xNewton; 
end

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
disp('%%%%%%%%%%%%%%%%%%PROBLEM #3 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
%% Part A
disp('%%%%%%%%Part 3A Solution:%%%%%%%');
% (a) Write a Python or MATLAB function that analytically solves a quadratic equation of the form:
% ax2 +bx+c = 0 
% and compare against roots that you find by hand via factorization or quadratic formula. 
%%use quadratic function to solve the column vector of coefficients [a, b, c] .
%%Test your code on the quadratic:
%%2x^2 − 6x+4=0 
syms x
coef=[2;-6;4];
f=@quadratic;
y=f(coef,x);
disp(y);
disp('%%%%%%%%End Part 3A Solution:%%%%%%%');
%% Part B, D, and E
disp('%%%%%%%%Part 3B Solution:%%%%%%%');
disp('Implementing Horners Method on Pg 194-195 to start, I was able expand it to calculate the derivative. ');
disp('However, I was not able to properly factor out the (x-5) term properly within the allotted time.');

% Part B Prompt:
% Write a polynomial division algorithm capable of dividing a given polynomial 
% Pn(x) (of order n and defined by a set of coefficients) by a given divisor (x − N).
% I.e. find Qn−1(x) such that: Pn(x) = (x − N)Qn−1(x) + R 
% Qn−1 is the polynomial left once you divide (x − N) out of Pn. Test your code by using it to
% divide out a factor of (x − 5) from the polynomial:
% x^5 −15x^4 +85x^3 −225x^2 +274x−120 = 0

%% Part C
disp('%%%%%%%%Part 3C Solution:%%%%%%%');
% Use an approximate Newton’s method (using an approximate derivative) to
% find one root of Eqn. 7: x^5 −15x^4 +85x^3 −225x^2 +274x−120 = 0
% A script to demonstrate solutions to nonlinear equations on closed
% and open domains
% requires:  objfun?.m (set function pointer f to desired function at beginning of program)
%% Params for Newton iteration
maxit=100;       %maximum number of iterations allowed
minx=0;
maxx=2*pi;
tol=1e-9;        %how close to zero we need to get to cease iterations

%% Objective function defs.
f=@objfun2;      %set the function for which we are finding roots, change to illustrate different problems
fprime=@objfun2_deriv;
x=linspace(minx,maxx,64);   %grid for basic plotting purposes
ygrid=f(x);
verbose=true;

%% Plot the function we are finding roots for
figure(2);
plot(x,ygrid);
title('Objective function')
xlabel('x')
ylabel('y')
axis tight;

%% Newton-Rhapson root-finding method
verbose=true;
[xNewton,itNew,flag]=newton_approx(f,-0.1*i,100,tol,verbose);
disp('1st Root value through Approx Newton method:  ');
disp(xNewton);
disp('Number of iterations required to reach tolerance:  ');
disp(itNew);

[xNewton,itNew,flag]=newton_approx(f,0.1*i,100,tol,verbose);
disp('1st Root value through Approx Newton method:  ');
disp(xNewton);
disp('Number of iterations required to reach tolerance:  ');
disp(itNew);

disp('%%%%%%%%End Part 3C Solution:%%%%%%%');

disp('%%%%%%%%Part 3D Solution:%%%%%%%');
% Part D Prompt: 
%Use your synthetic division algorithm to factor the root found in part c out of Eqn. 7 
%to produce a lower-order polynomial (i.e. Qn−1(x)). 
%The result will be a set of coefficients for a fourth order polynomial, denoted Q4.
disp('I was not able to calculate the result of the fourth order polynomial in the alotted time.');
disp('%%%%%%%%Part 3E Solution:%%%%%%%');
disp('I was able to properly "deflate" the polynomial to find the roots without using the quadratic solver.');
% Part E Prompt:
% Write a program to repeat the “deflation” process described above (parts c and d) 
% to find all of roots of Eqn. 7. This can be done by taking Qn-1 
% You can code your script to work specifically for this fifth order polynomial example
% it does not have to be general enough to work with a polynomial of arbitrary degree 
disp('The algorithm calculates the Five Roots of the polynomial. ');
soln1=@polySolver;
y = soln1(A,r0,nmax);
disp(y);
disp('%%%%%%%%End Part 3E Solution:%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%PROBLEM #3 ANSWER END%%%%%%%%%%%%%%%%%%');