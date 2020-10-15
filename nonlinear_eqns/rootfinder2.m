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

f=@objfunProblem2b;      
fprime=@objfunProblem2b_deriv;
verbose=true;
j=0;
finalarray3=[];

for g = -2:0.15:2 % changed variable name to g because of complex component
    [xNewton,itNew,flag]=newton_exactwithComplex(f,fprime,g.*i,maxit,tol,verbose);    
    j=j+1; 
    finalarray3(j)=xNewton; 
end
disp(finalarray3);
result1=finalarray3(1,6);
result2=finalarray3(1,26);
result3=finalarray3(1,20);
 
fprintf('The Real and Complex Roots of the polynomial are: %f, %f + %f i, and %f - %f i \n',real(result3),(result1),(result1),(result2),(result2));





