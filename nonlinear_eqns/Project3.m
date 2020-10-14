%% Author: Alec Sblendorio 
%% Due Date: October 15, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% Project #3
% Problem #1: Finding roots of functions lacking a closed form
% Problem #2: Numerical solution for multiple polynomial roots
% Problem #3: Multivariate root finding
%% Import test data for calibration of software

%% Problem #1: Finding roots of functions lacking a closed form
%%Part A:  Alter the Newton method function from the repository (newton exact.m)
%%so that it implements the approximate Newton method, i.e. so that the derivative 
%%is computed numerically as in Equation 3.77 in the course textbook.

%%Part B: Write a block of code that uses your newton approx.m function to 
%%find the first root (i.e. the smallest one) of the Bessel function of order zero:
%%J0(x) for the region 0 ≤ x ≤ ∞. Plot the Bessel function first and
%%use the plot to select a starting point in the vicinity of the first root.

%%Part C: Product a version of this script that finds the first six roots of the Bessel function
%%of order zero (these roots are needed for solutions of various types of 
%%boundary value problems in physics and engineering). 
%%These roots are commonly listed in ODE and PDE textbooks; 
%%look them up and verify your solutions (cite your sources).

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%Part A, B, C Solution:%%%%%%%');
soln = @Newton_Rhapson; 
x = soln();
disp('%%%%%%%%End Part A, B, C Solution:%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER END%%%%%%%%%%%%%%%%%%');
%% Problem #2: Numerical solution for multiple polynomial roots
disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
%%Part A: Suppose you have a polynomial of known or given order and need to find all of its roots.
%%Write a block of code or a script that uses the exact Newton method (e.g. the function in the course repository)
%%to find all of the real-valued roots of a polynomial. 
%%Assume that you do not need to identify repeated roots (if any). 
%%Test your code on the polynomial used in the book (Eqn. 3.115).

disp('%%%%%%%%Part A Solution:%%%%%%%');

disp('%%%%%%%%End Part A Solution:%%%%%%%');

%%Part B: Produce an altered version of your code to deal with the fact that there are potentially complex roots
%%to your polynomial. Use this code to find all roots, including complex-valued solutions
%%of the following polynomial (Eqn 3.145 in the book). 


disp('%%%%%%%%Part B Solution:%%%%%%%');

disp('%%%%%%%% End Part B Solution:%%%%%%%');


disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER END%%%%%%%%%%%%%%%%%%');

%% Problem #3: Multivariate root finding
disp('%%%%%%%%%%%%%%%%%%PROBLEM #3 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
%%Part A: Use the multi-dimensional Newton method
%%(given the in course repository: newton2D exact.m) to find all four roots of the system

disp('%%%%%%%%Part A Solution:%%%%%%%');
soln5 = @Newton_Rhapson2DAS;
s = soln5();

disp('%%%%%%%%End Part A Solution:%%%%%%%');

%%Part B: Produce an altered multi-dimensional Newton method (start from newton2D exact.m)
%%to find a root for the three equations system

disp('%%%%%%%%Part B Solution:%%%%%%%');

disp('%%%%%%%% End Part B Solution:%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #3 ANSWER END%%%%%%%%%%%%%%%%%%');
%% END PROJECT #3







