%% Author: Alec Sblendorio 
%% Due Date: October 8, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% Project #2
% Problem #1: LU factorization and its application to solve linear systems
% Problem #2: Iterative methods for solving linear systems

%% Import test data for calibration of software
Data1 = importdata('/Users/alecsblendorio/Documents/Projects/EP501_assignments/assignments/HW1/lowertriang_testproblem.mat');
Data2 = importdata('/Users/alecsblendorio/Documents/Projects/EP501_assignments/assignments/HW1/testproblem.mat');
%Data3 = importdata('/Users/alecsblendorio/Documents/Projects/EP501_assignments/assignments/HW2/iterative_testproblem.mat');
load('iterative_testproblem.mat');
A = Data2.A;
b = Data2.b;
b2 = Data2.b2;
b3 = Data2.b3;
L = Data1.L;
bL = Data1.bL;
% Ait = Data3.Ait;
% bit = Data3.bit;

%% Problem #1: LU factorization and its application to solve linear systems
%%Part A and B: New version of forward elimination to perform Doolittle-LU factorization
%%Using just the output of the factorization and a back-substitution function 
%%Solve the test linear system of equations given in the test problem
disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
soln = @DoolittleLU;
x = soln(A,b);
%%Part C: Use your LU factorized test matrix to set up a solution 
%%for the this system with different right hand sides


disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER END%%%%%%%%%%%%%%%%%%');
%% Problem #2: Iterative methods for solving linear systems
disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
%%Part A and B: Starting with the Jacobi function source code from the 
%repository create a new function that implements successive over-relaxation.
%Try this solver on the iterative test problem in this directory 
%Show that it gives the same results as the built-in Matlab utilities. 
% nit=10;
% Ait=diag(-1*ones(nit-1,1),-1)+diag(-1*ones(nit-1,1),1)+diag(4*ones(nit,1),0);    %this must be diagonally dominant or else the method won't converge
% %Ait=randn(nit,nit);    %see if code can detect non-diagonal dominance and exit gracefully...
% x0=randn(nit,1);
% bit=ones(nit,1);
% tol=1e-9;
% disp('Verbose Jacobi iterations:  ')
% %[xit,iterations]=sorFunc(x0,Ait,bit,tol,true);
% omega=1.10;
% [xit,iterations]=sorFunc(x0,Ait,bit,tol,omega,true);
% disp('%%%%%%%%Part A and B Solution:%%%%%%%');
% disp('Solution with Jacobi iterations:  ')
% disp(xit);
% disp('Number of iterations required and tolerance:  ')
% disp(iterations);
% disp(tol);
disp('%%%%%%%%Part A and B Solution:%%%%%%%');
disp('Solution with Jacobi iterations:  ')
disp(xit);
disp('Number of iterations required and tolerance:  ')
disp(iterations);
disp(tol);

disp('Matlab built-in solution:  ')
disp(Ait\bit);
disp('%%%%%%%%End Part A and B Solution:%%%%%%%');

%%Part C and D: By repeated application of your solver, iteratively adjust the relaxation 
%parameter until you find the approximate value that minimizes the number of iterations needed to achieve converge. 
%Use a fairly strict convergence criteria, e.g. 10−6 or less.
%How many fewer iterations are needed for your optimal case vs. t
%the standard Gauss-Seidel algorithm (relaxation parameter of one).
%Approximately speaking which values of relaxation parameter (if any) 
%perform worse than Gauss-Seidel?
disp('%%%%%%%%Part C Solution:%%%%%%%');

j=1;
tol=1e-9;
for i=0.25:0.05:1.5
    x=length(bit);
    x0=zeros(x,1);
    omega=i;
    [xit,iterations]=sorFunc(x0,Ait,bit,tol,omega,true);
    finalarray(1,j)=omega;
    finalarray(2,j)=iterations;
    j=j+1;
end %for 

disp('For Part C, the relaxation parameter is changed until the optimal value of iterations is reached.'); 
disp('Displayed in the matrix titled FinalArray, are various iterations that show the'); 
disp('optimal relaxation parameter.');
disp(finalarray); 
disp('%%%%%%%% End Part C Solution:%%%%%%%');

disp('%%%%%%%% Part D Solution:%%%%%%%');
disp('Number of iterations with Standard Gauss-Seidel Algorithm with a relaxtion parameter of 1:');
disp(24.0);
disp('Number of iterations with Optimal Case (relax parameter of 1.10)');
disp(19.0);

disp('Approximately, anything less than 1.0 and greater than 1.25 will perform worse than Gauss-Seidel.');

disp('%%%%%%%% End D Solution:%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER END%%%%%%%%%%%%%%%%%%');

%% END PROJECT #1 







