% A script to demonstrate solutions to nonlinear equations on closed
% and open domains
%
% requires:  objfun?.m (set function pointer f to desired function at beginning of program)

%% Params for Newton iteration
maxit=100;       %maximum number of iterations allowed
minx=0;
maxx=pi/4;
tol=1e-6;        %how close to zero we need to get to cease iterations

%% Objective function defs.
f=@objfun;      %set the function for which we are finding roots, change to illustrate different problems
fprime=@objfun_deriv;
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
% [xNewton,itNew,flag]=newton_exact(f,fprime,-1*i,100,tol,verbose);
% disp('Root value through Newton method:  ');
% disp(xNewton);
% disp('Number of iterations required to reach tolerance:  ');
% disp(itNew);
% 
% [xNewton,itNew,flag]=newton_exact(f,fprime,1*i,100,tol,verbose);
% disp('Root value through Newton method:  ');
% disp(xNewton);
% disp('Number of iterations required to reach tolerance:  ');
% disp(itNew);

j=0;
rec = 0;
finalarray1=[];

for i = 1:0.15:10
    [xNewton,itNew,flag]=newton_exact(f,fprime,10,100,tol,verbose);    
    j=j+1; 
    finalarray(j)=xNewton; 
end
 
disp(finalarray);
% result1=finalarray1(1,2);
% result2=finalarray1(1,3);
% result3=finalarray1(1,4);
% result4=finalarray1(1,5);
% result5=finalarray1(1,6);
% fprintf('The Five Roots of the polynomial are: %d,%d,%d,%d,%d\n',result1,result2,result3,result4,result5);
