%% Problem 2 (b)
%%Produce an altered version of your code to deal with the fact that there 
%%are potentially complex roots to your polynomial. Use this code to find all roots, 
%%including complex-valued solutions of the following polynomial
%%x3 −3x2 +4x−2 = 0

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

 result1=finalarray3(1,1);
 result2=finalarray3(1,11);
 
 fprintf('The Real and Complex Roots of the polynomial are: %f + %f i and %f - %f i \n',real(result1),imag(result1),real(result2),imag(result2));





