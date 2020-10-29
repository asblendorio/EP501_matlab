% Problem 1a 
% Write a program that performs a linear least squares fit of a set of data to a polynomial
% of arbitrary order n. You may use any functions in the course repository 
% or that you have written for your homework for solving the required system of equations
% however, you may not use the built-in MATLAB functions for your solution (only to check the results).
% Problem 1b 
% Use your fitting program to fit the test data to a line and a quadratic form.
% Plot your results and the data on the same axis so they can be easily compared. 
% Test your results against the built-in Matlab functions polyfit and polyval 
% Compare the error vectors and residuals for the two fits.

function [x] = lesq(ndim,ndeg,nref,x,f,a,b,c)
nref = length(); % number of data points and polynomial coefficients
x = size(10,1);
f = size(10,1);
a = size(10,10);
b = size(10,1);
c = size(10,1); 
sumx=0.0;
sumxx=0.0;
sumf=0.0;
sumxf=0.0;
for i = 1:nref
    sumx=sumx+x(i);
    sumxx=sumxx+x(i).*2;
    sumf=sumf+f(i);
    sumxf=sumxf+x(i).*f(i);
end %for

a(1,1) = nref;
a(1,2) = sumx;
a(2,1) = sumx;
a(2,2) = sumxx;
b(1) = sumf;
b(2) = sumxf;

%Perform Gaussian Elimination here 
verbose = true;
gaussx = Gauss_elim(A,b,true);
disp(gaussx);

%% Plotter 


end %function 