%%Problem #3a
%%Write a Python or MATLAB function that analytically solves a quadratic
%%equation of the form: %  A(x^2)+B(x)+C=0
%%and compare against roots that you find by hand via factorization or quadratic formula. 
function [s,j]=quadratic(coef,x)
q = input('Press 1 to show general solution. Press 2 to solve analytically.');

if q == 1
syms a b c x
eqn1 = a.*x.^2 + b.*x + c == 0;
s = solve(eqn1);
disp('The solution to the quadratic equation is:');
disp(s);
end %if

if q == 2
eqn2 = (coef(1,1).*x.^2) +(coef(2,1).*x)+(coef(3,1)) == 0;
j = solve(eqn2);
disp('The roots to the input quadratic equation is:');
disp(j);
end %if 

end %function
