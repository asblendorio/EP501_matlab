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

function [coeffs] = lesq(x,ynoisy,sigmay,n)
nref = length(x); % number of data points and polynomial coefficients
J=cat(2,ones(nref,1),x(:)); 
if n >= 1
    for i=2:n
        J=cat(2,J,x(:).^i);
    end %for
end %if

M=J'*J;
yprime=J'*ynoisy(:);
[Mmod,ord]=Gauss_elim(M,yprime);
avec=backsub(Mmod(ord,:));
% Manual Linear Fit 
yfit=avec(1);
% yfitlin=avec(1)+avec(2).*x;
% yfitquad=avec(1);
% loop over yfit data for quadratic and higher order 
for i=1:n
    yfit=yfit+avec(i+1).*x.^i;
end %for

% calculate the error and residual from Dr. Z
% er = abs(yfit-ynoisy);  
% residual = sum(er);
% disp('Error');
% disp(residual);

%% Plotter and Compare against MATLAB and Polyval
figure(1);
plot(x,ynoisy,'o','MarkerSize',2);
hold on; 

l=n;
coeffs=polyfit(x,ynoisy,l);
xlarge=linspace(-1,1,50);
ylarge=polyval(coeffs,xlarge);
hold on;

plot(xlarge,ylarge,'k--','linewidth',2);
hold on;
%linear fit 

plot(x,yfit,'r--','linewidth',2);
legend('Input Data','True Matlab Function ','New Algorithm Fit');
hold on;

fprintf('Coefficients of polynomial of order ,%d, from Linear Least Square Algorithm:',(n));
disp(avec);
disp('Matlab,GNU/Octave built-in solution:');
disp(coeffs);

%% Perform Chi Squared Goodness of Fit 
[chi2] = chi_squared(ynoisy,yfit,n,sigmay);
disp(chi2);
end %function 