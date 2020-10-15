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
f=@objfunProblem2;      %set the function for which we are finding roots, change to illustrate different problems
fprime=@objfunProblem2_deriv;
verbose=true;
j=0;
finalarray2=[];

for i = 0:1.0:10
    [xNewton,itNew,flag]=newton_exact(f,fprime,i,maxit,tol,verbose);    
    j=j+1; 
    finalarray2(j)=xNewton; 
end

result1=finalarray2(1,2);
result2=finalarray2(1,3);
result3=finalarray2(1,4);
result4=finalarray2(1,5);
result5=finalarray2(1,6);
fprintf('The Roots of the polynomial are: %d,%d,%d,%d,%d',result1,result2,result3,result4,result5);





