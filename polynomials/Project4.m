%% Author: Alec Sblendorio 
%% Due Date: November 6, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% Project #4
%%Problem 1: This problem concerns least squares and data fitting and requires use of the example dataset from the repository, test lsq.mat, which provides data for variables xi, yi, σyi referenced below.
%%Problem #2: This problem concerns bilinear interpolation methods and requires use of the grid (variables xg,yg) and data samples (f2D) from test interp.mat
%% Data Input 
load('test_lsq.mat');
load('test_interp.mat');
%% Problem #1: Finding roots of functions lacking a closed form
%%Part A: Write a program that performs a linear least squares fit of a set 
%%of data to a polynomial of arbitrary order n. You may use any functions 
%%in the course repository or that you have written for your homework for 
%%solving the required system of equations; however, you may not use the 
%%built-in MATLAB functions for your solution (only to check the results).

%%Part B: Use your fitting program to fit the test data (yi located at independent variable positions xi) 
%%to a line and a quadratic form. Plot your results and the data on the 
%%same axis so they can be easily compared. Test your results against the
%%built-in Matlab functions polyfit and polyval. Compare the error vectors 
%%and residuals for the two fits.

%%Part C: A rigorous way of deciding between preferred functional forms in fits
%%(from some set of options like linear, quadratic, cubic, quartic, etc.)
%%is to define a goodness-of-fit statistic that quantifies how effective a 
%%particular form is a describing a given data set. This goodness statistic 
%%should balance the need to fit the data (by having more parameters (unknowns)
%%in the fit against the fact that of course one can fit a set of data given
%%enough unknowns (a problem referred to as “overfitting the data”). 
%%The simplest and most commonly used goodness of fit statistic is the
%%reduced Chi-squared statistic. Write a function that evaluates χ2ν for a 
%%polynomial fit of order n.

%%Part D: Use your goodness-of-fit statistic to determine whether is best 
%%to fit these data with a linear, quadratic, or cubic polynomial.
%%Show how you reached your decision.

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
disp('This is for linear fit');
soln1=@lesq;
a=soln1(x,ynoisy,sigmay,1);
disp(a);

disp('This is for quadratic fit');
soln2=@lesq;
b=soln2(x,ynoisy,sigmay,2);
disp(b);

disp('This is for cubic fit');
soln3=@lesq;
c=soln3(x,ynoisy,sigmay,3);
disp(c);

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER END%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER PLOTS%%%%%%%%%%%%%%%%%%'
%% Problem #2: This problem concerns bilinear interpolation methods and requires use of the grid (variables xg,yg) and data samples (f2D) from test interp.mat.
%%Part A:  Write a function that takes in a grid of points describing some independent variable
%%(say xi), and a point to which the data are to be interpolated x′ and 
%%finds the index i into the array xi such that: xi ≤ x′ ≤ xi+1.

%%Part B: Use the function from part (a) to construct an additional function 
%%that works over a 2D grid x, y. I.e. given two grids xi, yj find the indices i, j 
%%such that: xi ≤ x′ ≤ xi+1, yj ≤ y′ ≤ yj+1.

%%Part C: Use your results from parts a and b to create a bilinear interpolation function
%%that takes in a sequence of data points {x′k,yk′ } to which data are being 
%%interpolated, a grid xi,yj, and a dataset fij that is defined over this grid
%%and produces bilinearly interpolated values of fk at the points {x′k,yk′ }.
%%Write your program so that the input points are simply a flat list and not
%%necessarily a 2D grid of points (you can always reshape the results later if needed).

%%Part D: Test your results against Matlab’s bilinear interpolation function (interp2)
%%and show that you get the same result. Use the test data from the repository 
%%(test interp.mat). The source grid data are stored in xg,yg, while the value 
%%of the function at those points is in f2D. xgi,ygi are the densely sampled 
%%grid point to which the data are to be interpolated for this test.
disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');

disp('%%%%%%%%Part 2A Solution:%%%%%%%');


disp('%%%%%%%%End Part 2A Solution:%%%%%%%');

disp('%%%%%%%%Part 2B Solution:%%%%%%%');

disp('%%%%%%%% End Part 2B Solution:%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER END%%%%%%%%%%%%%%%%%%');


%% END PROJECT #4