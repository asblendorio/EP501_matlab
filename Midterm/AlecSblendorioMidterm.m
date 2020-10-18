%% Author: Alec Sblendorio 
%% Due Date: October 19, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% MIDTERM EXAM 
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


%% Problem 2

%% Problem 3 
%% Part A
% (a) Write a Python or MATLAB function that analytically solves a quadratic equation of the form:
% ax2 +bx+c = 0 
% and compare against roots that you find by hand via factorization or quadratic formula. 
%%use quadratic function to solve the column vector of coefficients [a, b, c] .
%%Test your code on the quadratic:
%%2x^2 − 6x+4=0 
syms x
coef=[2;-6;4];
f=@quadratic;
x=f(coef,x);
disp(x);
%% Part B
% Write a polynomial division algorithm capable of dividing a given polynomial 
% Pn(x) (of order n and defined by a set of coefficients) by a given divisor (x − N).
% I.e. find Qn−1(x) such that: Pn(x) = (x − N)Qn−1(x) + R 
% Qn−1 is the polynomial left once you divide (x − N) out of Pn. Test your code by using it to
% divide out a factor of (x − 5) from the polynomial:
% x^5 −15x^4 +85x^3 −225x^2 +274x−120 = 0

%% Part C
% Use an approximate Newton’s method (using an approximate derivative) to
% find one root of Eqn. 7
% x^5 −15x^4 +85x^3 −225x^2 +274x−120 = 0
% A script to demonstrate solutions to nonlinear equations on closed
% and open domains
%
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

%% Plot the function we are finding roots for
figure(1);
plot(x,ygrid);
title('Objective function')
xlabel('x')
ylabel('y')
axis tight;


%% Newton-Rhapson root-finding method
verbose=true;
j=0;
rec = 0;
finalarray1=[];

for i = 0:1.0:10
    [xNewton,derivative]=newton_approx(f,i,maxit,tol,verbose);    
    j=j+1; 
    finalarray1(j)=xNewton; 
end

result1=finalarray1(1,2);
result2=finalarray1(1,3);
result3=finalarray1(1,4);
result4=finalarray1(1,5);
result5=finalarray1(1,6);
fprintf('The Five Roots of the polynomial are: %d,%d,%d,%d,%d\n',result1,result2,result3,result4,result5);

%% Part D 







