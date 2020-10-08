%% Author: Alec Sblendorio 
%% Due Date: October 8, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% Project #2
% Problem #1: LU factorization and its application to solve linear systems
% Problem #2: Iterative methods for solving linear systems

%% Import test data for calibration of software
Data1 = importdata('/Users/alecsblendorio/Documents/Projects/EP501_assignments/assignments/HW1/lowertriang_testproblem.mat');
Data2 = importdata('/Users/alecsblendorio/Documents/Projects/EP501_assignments/assignments/HW1/testproblem.mat');

A = Data2.A;
b = Data2.b;
b2 = Data2.b2;
b3 = Data2.b3;
L = Data1.L;
bL = Data1.bL;

%% Problem #1: LU factorization and its application to solve linear systems
%%Part A and B: New version of forward elimination to perform Doolittle-LU factorization
disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
soln = @DoolittleLU;
x = soln(A,b);


disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER END%%%%%%%%%%%%%%%%%%');
%% Problem #2: Iterative methods for solving linear systems
disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
%%Part A:   
%%B,C, and D: 
soln3 = @gaussjordanElim;
g = soln3(A,b);
disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER END%%%%%%%%%%%%%%%%%%');

%% END PROJECT #1 







