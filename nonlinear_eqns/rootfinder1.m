%% Problem 2 (a)
%%Suppose you have a polynomial of known or given order and need to find all of its roots.
%%Write a block of code or a script that uses the exact Newton method 
%% Params for Newton iteration
maxit=100;       %maximum number of iterations allowed
minx=0;
maxx=2*pi;
tol=1e-9;        %how close to zero we need to get to cease iterations

%% Objective function defs.
f=@objfun2;      %set the function for which we are finding roots, change to illustrate different problems
fprime=@objfun2_deriv;
x=linspace(minx,maxx,64);   %grid for basic plotting purposes
ygrid=f(x);
verbose=true;

j=0;
rec = 0;
finalarray2=[];

for i = 0:0.5:10
    [xNewton,itNew,flag]=newton_exact(f,fprime,i,maxit,tol,verbose);  % The only difference between repository's and mine is the iteration records not outputted. 
    j=j+1;  
    finalarray2(j)=xNewton;       % this is a data bank of the roots computed from the iteration process, may have repeats of the same root, but does not matter -- just pick unique ones for the final results. 
end

sorted=sort(finalarray2);
disp('between my guess of 0 to 10 with increment of 0.5, the roots I could obtain are: ')
disp(sorted);

disp('Therefore, the roots of the equation 3.115, without identity repeated roots, are: ')
disp('1, 2, 3, 5, and 5')



