% A script to demonstrate solutions to nonlinear equations on closed
% and open domains
%
% requires:  objfun?.m (set function pointer f to desired function at beginning of program)


%% Params for Newton iteration
maxit=100;       %maximum number of iterations allowed
minx=0;
maxx=pi/4;
tol=1e-9;        %how close to zero we need to get to cease iterations

%% Objective function defs.
f=@objfun;      %set the function for which we are finding roots, change to illustrate different problems
fprime=@objfun_deriv;
y = polyval(p,x);
x=linspace(minx,maxx,64);   %grid for basic plotting purposes
ygrid=f(x);
verbose=true;

%% Plot the function we are finding roots for
figure(1);
plot(x,ygrid);
title('Objective function')
xlabel('x')
ylabel('y')
axis tight;


%% Newton-Rhapson root-finding method
verbose=true;
j=0;
rec = 0;
finalarray1=[];

for i = 0:1.0:10
    [xNewton,derivative]=newton_approx(f,i,maxit,tol,verbose);    
    j=j+1; 
    finalarray1(j)=xNewton; 
end

result1=finalarray1(1,2);
result2=finalarray1(1,3);
result3=finalarray1(1,4);
result4=finalarray1(1,5);
result5=finalarray1(1,6);
fprintf('The Five Roots of the polynomial are: %d,%d,%d,%d,%d\n',result1,result2,result3,result4,result5);
