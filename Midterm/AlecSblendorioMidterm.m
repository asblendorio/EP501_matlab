%% Author: Alec Sblendorio 
%% Due Date: October 19, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% MIDTERM EXAM 
% Problem #1: In class we compared the efficiency of Gaussian elimination vs. Jacobi iteration for linear systems of various size in the benchmark scripts (see ./linear algebra/benchmark.(py,m) in the course code repository for the language of your choice). In this problem we will explore a special method for solving tridiagonal systems that is even more efficient than these two options - the Thomas algorithm.
% (a) Write a MATLAB or Python function tridiag() that solves a tridiagonal system of equations using the Thomas algorithm.
% Verify your solution by applying it to the iterative test problem for HW2 in the EP501 Assignments repository 
% in ./HW2/iterative testproblem.mat file (this problem is actually tridiagonal). 
% (b) Produce a new version of the benchmark.(py,m) script that runs timed tests for all three of: Gaussian elimination, Jacobi iteration, and your new tridiag solver. Use the same tridiagonal problem for testing each of the solvers (as in the existing benchmark.(py,m) script). Store the solve times for each method and each system size in separate arrays and plot times vs. number of unknowns for each method on a single graph.
% (c) For what system sizes does Gaussian elimination perform better than Jacobi iteration (on your computer)?
% (d) Explain why the Thomas algorithm is more efficient for large numbers of unknowns. Try to be as quantitative as possible by discussing the number of operations requires to execute each algorithm/iteration. 