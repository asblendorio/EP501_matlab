%% Author: Alec Sblendorio 
%% Due Date: October 21, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% MIDTERM EXAM 
%% Import Data 
load('iterative_testproblem.mat');
load('PolynomialTest.mat');
%% Problem 1 
% Problem #1: In class we compared the efficiency of Gaussian elimination vs. Jacobi iteration for linear systems of various size in the benchmark scripts 
% In this problem we will explore a special method for solving tridiagonal systems that is even more efficient than these two options - the Thomas algorithm.
% (a) Write a MATLAB or Python function tridiag() that solves a tridiagonal system of equations using the Thomas algorithm.
% Verify your solution by applying it to the iterative test problem for HW2 in the EP501 Assignments repository 
% in ./HW2/iterative testproblem.mat file (this problem is actually tridiagonal). 
% (b) Produce a new version of the benchmark.(py,m) script that runs timed tests for all three of: 
% Gaussian elimination, Jacobi iteration, and your new tridiag solver. Use the same tridiagonal problem for 
% testing each of the solvers (as in the existing benchmark.(py,m) script). Store the solve times for each method and each system size in separate arrays and plot times vs. number of unknowns for each method on a single graph.
% (c) For what system sizes does Gaussian elimination perform better than Jacobi iteration (on your computer)?
% (d) Explain why the Thomas algorithm is more efficient for large numbers of unknowns. 
% Try to be as quantitative as possible by discussing the number of operations requires to execute each algorithm/iteration.
%% Part A 
% (a) Write a MATLAB or Python function tridiag() that solves a tridiagonal system of equations using the Thomas algorithm.
% Verify your solution by applying it to the iterative test problem for HW2 in the EP501 Assignments repository 
% in ./HW2/iterative testproblem.mat file (this problem is actually tridiagonal).
%Data = importdata('/Users/alecsblendorio/Documents/Projects/EP501_matlab/Midterm/iterative_testproblem.mat');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%Part 1A Solution:%%%%%%%');
soln = @tridiag; 
x = soln(Ait,bit);
disp('Thomas Algorithm Solution:');
disp(x);
disp('%%%%%%%%End Part 1A Solution:%%%%%%%');
%% Part B
%%Produce a new version of the benchmark.(py,m) script that runs timed tests for all three of:
%%Gaussian elimination, Jacobi iteration, and your new tridiag solver. Use the same tridiagonal
%%problem for testing each of the solvers (as in the existing benchmark.(py,m) script). Store the
%%solve times for each method and each system size in separate arrays and plot times vs. number
%%of unknowns for each method on a single graph.

% Evaluate performance and scaling of Gaussian elimination, Jacobi iteration,
% and Tri-Diagonal Solver by solving systems of different size and timing the solves
disp('%%%%%%%%Part 1B Solution:%%%%%%%');
nvals=50:50:500;
testtimes=zeros(size(nvals));
lrep=10;     %how many times to repeat each test

disp('Start of tests of Gaussian-elimination scaling');
for in=1:numel(nvals)
    nlarge=nvals(in);
    Blarge=diag(-1*ones(nlarge-1,1),-1)+diag(-1*ones(nlarge-1,1),1)+diag(4*ones(nlarge,1),0);    %this must be diagonally dominant or else the method won't converge
    blarge=ones(nlarge,1);
    
    for irep=1:lrep     %benchmark will repeat the same solution several times to eliminate random variations from CPU load, etc.
        tstart=cputime;
        [Blargemod,ordlarge]=Gauss_elim(Blarge,blarge);
        xlarge=backsub(Blargemod(ordlarge,:));
        tend=cputime;
        testtimes(in)=testtimes(in)+(tend-tstart)/lrep;
    end %for
    disp([' GE solution for system of size ',num2str(nlarge),' takes average time ',num2str(testtimes(in)),' s']);
end %for

figure(1);
plot(nvals,testtimes,'o','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','blue')
xlabel('system size');
ylabel('time to solve (s)');
title('Empirically Determined Performance');

disp('Start of tests for Jacobi iteration');
tol=1e-9;
testtimes=zeros(size(nvals));
for in=1:numel(nvals)
    nlarge=nvals(in);
    Blarge=diag(-1*ones(nlarge-1,1),-1)+diag(-1*ones(nlarge-1,1),1)+diag(4*ones(nlarge,1),0);    %this must be diagonally dominant or else the method won't converge
    blarge=ones(nlarge,1);

    for irep=1:lrep     %benchmark will repeat the same solution several times to eliminate random variations from CPU load, etc.
        tstart=cputime;
        x0=randn(nlarge,1);
        [xit,iterations]=Jacobi(x0,Blarge,blarge,tol,false);
        tend=cputime;
        testtimes(in)=testtimes(in)+(tend-tstart)/lrep;
    end %for
    disp([' JI solution for system of size ',num2str(nlarge),' takes average time ',num2str(testtimes(in)),' s']);
end %for

figure(1);
hold on
plot(nvals,testtimes,'^','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','blue')
xlabel('system size');
ylabel('time to solve (s)');
legend('Gauss elim.','Jacobi it.')
title('Empirically Determined Performance');

disp('Start of tests for Tridiagonal');
for in=1:numel(nvals)
    nlarge=nvals(in);
    Blarge=diag(-1*ones(nlarge-1,1),-1)+diag(-1*ones(nlarge-1,1),1)+diag(4*ones(nlarge,1),0);    %this must be diagonally dominant or else the method won't converge
    blarge=ones(nlarge,1);
    for irep=1:lrep     %benchmark will repeat the same solution several times to eliminate random variations from CPU load, etc.
        tstart=cputime;
        [Awork]=tridiag(Blarge,blarge);
        xlarge=backsub(Awork(1,:));
        tend=cputime;
        testtimes(in)=testtimes(in)+(tend-tstart)/lrep;
    end %for
    disp([' Tridiagonal solution for system of size ',num2str(nlarge),' takes average time ',num2str(testtimes(in)),' s']);
    
end %for

figure(1);
hold on
plot(nvals,testtimes,'*','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','blue')
xlabel('system size');
ylabel('time to solve (s)');
legend('Gauss elim.','Jacobi it.','TriDiag')
title('Empirically Determined Performance');
disp('%%%%%%%%End Part 1B Solution:%%%%%%%');
%% Part C & D
% % Handwritten Components are scanned and added in here
disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER END%%%%%%%%%%%%%%%%%%');
%% Problem 2
disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
% Handwritten Components are scanned and added in here
%% Part A 
% Under what condition does this homogeneous system of equations 
% have a nontrivial solution (viz. u = 0)?
%% Part B 
% Derive a polynomial for the unknown wave speed v using the condition you found in part a of this problem.
%% Part C 
% How many roots (known as wave modes) does your characteristic polynomial have?
%% Part D 
% Two of the roots for this system are v = ±Ca cos θ. 
% Plug one of these roots back into the matrix equation and use it to determine
% which component(s) of drift, i.e. ux,y,z, can be nonzero for this wave mode. 
%% Part E 
% Given typical parameters for the solar wind plasma:
% γ = 5/3, ρ = 1.67 × 10−21 [kg/m3], p = 1.38 × 10−11 [Pa], and B = 10−9 [T]
% compute values for the sound and Alfven speeds in the plasma.
%% Part F
disp('%%%%%%%%Part 2F Solution:%%%%%%%');
disp('Part 2F has required considerable time to figure out a sensible convergence.');
disp('I was not able to identify the correct criteria within the timeframe.');
% Use an exact Newton method (i.e. use analytically computed derivatives)
% to find numerical values for all roots using θ = π/4 for the angle of propagation. 
% Make sure you select a sensible convergence criteria given the coefficients 
% for this problem and treat all parameters except for the unknown roots for v 
% to be constant for purposes of developing derivatives needed to implement Newton’s method.

% A script to demonstrate solutions to nonlinear equations on closed
% and open domains
% requires:  objfun?.m (set function pointer f to desired function at beginning of program)

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
disp('%%%%%%%%Part 3B, D, and E Solution:%%%%%%%');
disp('Implementing Horners Method on Pg 194-195 to start, I was able expand it to calculate the derivative. ');
disp('The algorithm calculates the Five Roots of the polynomial. ');
% Part B Prompt:
% Write a polynomial division algorithm capable of dividing a given polynomial 
% Pn(x) (of order n and defined by a set of coefficients) by a given divisor (x − N).
% I.e. find Qn−1(x) such that: Pn(x) = (x − N)Qn−1(x) + R 
% Qn−1 is the polynomial left once you divide (x − N) out of Pn. Test your code by using it to
% divide out a factor of (x − 5) from the polynomial:
% x^5 −15x^4 +85x^3 −225x^2 +274x−120 = 0

% Part D Prompt: 
%Use your synthetic division algorithm to factor the root found in part c out of Eqn. 7 
%to produce a lower-order polynomial (i.e. Qn−1(x)). 
%The result will be a set of coefficients for a fourth order polynomial, denoted Q4.

% Part E Prompt:
% Write a program to repeat the “deflation” process described above (parts c and d) 
% to find all of roots of Eqn. 7. This can be done by taking Qn-1 
% You can code your script to work specifically for this fifth order polynomial example
% it does not have to be general enough to work with a polynomial of arbitrary degree 

soln1=@polySolver;
y = soln1(A,r0,nmax);
disp(y);

disp('%%%%%%%%Part 3B, D, and E Solution:%%%%%%%');
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

disp('%%%%%%%%%%%%%%%%%%PROBLEM #3 ANSWER END%%%%%%%%%%%%%%%%%%');