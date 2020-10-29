% Problem 1c 
% Write a function that evaluates reduced chi^2 value for a polynomial fit of order n.
% x_i is the points of the independent variable at which data are sampled,
% y_i is the data
% f(x) is the function which is being fitted to the data (linear, quadratic, cubic, etc.)
% σ_i is the uncertainty of the data 
% ν is the number of degrees of freedom in your fit.
% ν  = number of data points(n) minus number of parameters(P) being estimated
function chi2 = chi_squared(ynoisy,x,P,sigmay)

n = max(size(ynoisy));
terms = (((ynoisy-x).^2)./sigmay).^2;
dof = (n-P);
chi2 = (1./dof)*sum(terms);
disp('degrees of freedom:');
disp(dof);
t = chi2/dof;
disp('testing if fit is good:');
disp(t);

disp('The reduced Chi Squared Value of a polynomial of order N is: ');

end %function 

