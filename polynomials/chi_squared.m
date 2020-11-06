% Problem 1c 
% Write a function that evaluates reduced chi^2 value for a polynomial fit of order n.
% x_i is the points of the independent variable at which data are sampled,
% y_i is the data
% f(x) is the function which is being fitted to the data (linear, quadratic, cubic, etc.)
% σ_i is the uncertainty of the data 
% ν is the number of degrees of freedom in your fit.
% ν  = number of data points(n) minus number of parameters(P) being estimated
function chi2 = chi_squared(ynoisy,yfit,P,sigmay)
n = max(size(ynoisy));
finalarray=zeros(length(n));
t=1;
for i=1:n
    terms = (((ynoisy(i)-yfit(i))).^2/(sigmay(i).^2));
    finalarray(1,t)=terms;
    t=t+1;
end %for
v = (n-P);
chi2 = (1./v).*sum(finalarray);

disp('The reduced Chi Squared Value of a polynomial of order:');

end %function 

