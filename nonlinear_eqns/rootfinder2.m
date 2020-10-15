%% Problem 2 (a)
%%Suppose you have a polynomial of known or given order and need to find all of its roots.
%%Write a block of code or a script that uses the exact Newton method 
%% Params for Newton iteration
maxit=100;       %maximum number of iterations allowed
minx=0;
maxx=2*pi;
tol=1e-9;        %how close to zero we need to get to cease iterations

%% Newton-Rhapson root-finding method for polynomials 
%%this method will find all of the REAL-valued roots of a polynomial
f=@objfunProblem2b;      %set the function for which we are finding roots, change to illustrate different problems
fprime=@objfunProblem2b_deriv;
x=linspace(minx,maxx,64);   %grid for basic plotting purposes
ygrid=f(x);
verbose=true;
j=0;
rec = 0;
finalarray3=[];

for g = -5:0.25:5 % changed variable name to g because of complex component
    [xNewton,itNew,flag]=newton_exact(f,fprime,g+g.*i,maxit,tol,verbose);    
    j=j+1; 
    finalarray3(j)=xNewton; 
end
disp(finalarray3);
% result1=finalarray3(1,2);
% result2=finalarray3(1,3);
% result3=finalarray3(1,4);
% result4=finalarray3(1,5);
% result5=finalarray3(1,6);
% fprintf('The Real and Complex Roots of the polynomial are: %d,%d,%d,%d,%d',result1,result2,result3,result4,result5);





