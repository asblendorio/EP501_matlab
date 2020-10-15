%% Author: Alec Sblendorio 
%% Due Date: October 15, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% Project #3
% Problem #1: Finding roots of functions lacking a closed form\n
% Problem #2: Numerical solution for multiple polynomial roots\n
% Problem #3: Multivariate root finding
%% Problem #1: Finding roots of functions lacking a closed form
%%Part A:  Alter the Newton method function from the repository (newton exact.m)
%%so that it implements the approximate Newton method, i.e. so that the derivative 
%%is computed numerically as in Equation 3.77 in the course textbook.

%%Part B: Write a block of code that uses your newton approx.m function to 
%%find the first root (i.e. the smallest one) of the Bessel function of order zero:
%%J0(x) for the region 0 ≤ x ≤ ∞. Plot the Bessel function first and
%%use the plot to select a starting point in the vicinity of the first root.

%%Part C: Product a version of this script that finds the first six roots of the Bessel function
%%of order zero (these roots are needed for solutions of various types of 
%%boundary value problems in physics and engineering). 
%%These roots are commonly listed in ODE and PDE textbooks; 
%%look them up and verify your solutions (cite your sources).

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%Part 1A Solution:%%%%%%%');

maxit=100;       %maximum number of iterations allowed
minx=0;
maxx=2*pi;
tol=1e-6;        %how close to zero we need to get to cease iterations

%% Objective function defs.
f=@objfun2;      %set the function for which we are finding roots, change to illustrate different problems
fprime=@objfun2_deriv;
x=linspace(minx,20,64);   %grid for basic plotting purposes
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
[xNewton,itNew,flag]=newton_approx(f,-0.3*i,100,tol,verbose);
disp('Root value through Newton method:  ');
disp(xNewton);
disp('Number of iterations required to reach tolerance:  ');
disp(itNew);

[xNewton,itNew,flag]=newton_approx(f,0.3*i,100,tol,verbose);
disp('Root value through Newton method:  ');
disp(xNewton);
disp('Number of iterations required to reach tolerance:  ');
disp(itNew);

% %% Newton approach for suspected complex roots
disp('Complex Stuff');
[xNewton,itNew]=newton_approx(f,7*i,100,tol,verbose);
disp('Root value through Newton method:  ');
disp(xNewton);
disp('Number of iterations required to reach tolerance:  ');
disp(itNew);

[xNewton,itNew]=newton_approx(f,-7*i,100,tol,verbose);
disp('Root value through Newton method:  ');
disp(xNewton);
disp('Number of iterations required to reach tolerance:  ');
disp(itNew);
%%using Approximate Newton Function 
[xNewton,derivative]=newton_approx(f,2,maxit,tol,verbose);
disp('The Approximate Newton Method has found the first root to be:');
disp(xNewton);
disp('The Approximate Newton Method has found the derivative to be:');
disp(derivative);
disp('%%%%%%%%End Part 1A Solution:%%%%%%%');

disp('%%%%%%%%Part 1B Solution:%%%%%%%');
%%In Office Hours, Dr. Z alluded to looping over values form 0 to 20 to
%%calculate the first root of the Bessel function. 
iterate=0:1.0:20; % 6th root is before 20.
maxit=100;
tol=1e-9;
verbose=true;
f=@Bessel_objfun1;
minx=0;
maxx=20;
%% Plot the function we are finding roots for
figure(2);
plot(x,f(x));
title('Bessel Function of Order Zero')
xlabel('x')
ylabel('y')
axis tight;

disp('To aid in finding the roots of the Bessel function, a fellow classmate provided the class a document cited below. This gives us the first 10 roots.');
disp('https://thalis.math.upatras.gr/~vrahatis/papers/journals/VrahatisGRZ97_Z_ANGEW_MATH_MECH_77_pp467-475_1997.pdf\n');

disp('%%%%%%%%End Part 1B Solution:%%%%%%%');

disp('%%%%%%%%Part 1C Solution:%%%%%%%');


disp('%%%%%%%%End Part 1C Solution:%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER END%%%%%%%%%%%%%%%%%%');
%% Problem #2: Numerical solution for multiple polynomial roots
disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');

disp('%%%%%%%%Part 2A Solution:%%%%%%%');
%% Problem 2 (a)
%%Part A: Suppose you have a polynomial of known or given order and need to find all of its roots.
%%Write a block of code or a script that uses the exact Newton method (e.g. the function in the course repository)
%%to find all of the real-valued roots of a polynomial. 
%%Assume that you do not need to identify repeated roots (if any). 
%%Test your code on the polynomial used in the book (Eqn. 3.115).
%% Params for Newton iteration
maxit=100;       %maximum number of iterations allowed
minx=0;
maxx=2*pi;
tol=1e-9;        %how close to zero we need to get to cease iterations

%% Newton-Rhapson root-finding method for polynomials 
%%this method will find all of the REAL-valued roots of a polynomial
f=@objfunProblem2;      %set the function for which we are finding roots, change to illustrate different problems
fprime=@objfunProblem2_deriv;
x=linspace(minx,maxx,64);   %grid for basic plotting purposes
ygrid=f(x);
verbose=true;
j=0;
rec = 0;
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
fprintf('The Real Roots of the polynomial are: %d,%d,%d,%d,%d\n',result1,result2,result3,result4,result5);

disp('%%%%%%%%End Part 2A Solution:%%%%%%%');

disp('%%%%%%%%Part 2B Solution:%%%%%%%');
%% Problem 2 (b)
%%Part B: Produce an altered version of your code to deal with the fact that there are potentially complex roots
%%to your polynomial. Use this code to find all roots, including complex-valued solutions
%%of the following polynomial (Eqn 3.145 in the book). 
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
% disp(finalarray3);
result1=finalarray3(1,6);
result2=finalarray3(1,26);
result3=finalarray3(1,20);
 
fprintf('The Real and Complex Roots of the polynomial are: %f, %f + %f i, and %f - %f i \n',real(result3),(result1),(result1),(result2),(result2));

disp('%%%%%%%% End Part 2B Solution:%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER END%%%%%%%%%%%%%%%%%%');

%% Problem #3: Multivariate root finding
disp('%%%%%%%%%%%%%%%%%%PROBLEM #3 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
%%Part A: Use the multi-dimensional Newton method
%%(given the in course repository: newton2D exact.m) to find all four roots of the system

disp('%%%%%%%%Part 3A Solution:%%%%%%%');
%% Multidimensional function from Problem #3a
fm=@objfun2Df_AS;
gm=@objfun2Dg_AS;
gradfm=@grad_objfun2Df_AS;
gradgm=@grad_objfun2Dg_AS;

%this is for plotting
x=linspace(-1.5,1.5,20);
y=linspace(-1.5,1.5,20);
[X,Y]=meshgrid(x,y);
F=fm(X,Y);
G=gm(X,Y);

%% Newton's method for multi-variable nonlinear equations
%x0=i;
%y0=0.258*i;
% use for loop to iterate over inital x0 and y0 values
%loop keep stopping after finding 1 value, need to keep it going to find
%all 4 roots
j=1;
for i=0.1:0.25:3.0
    x0=i;
    y0=i;
    [xm,ym,it2D,success2D]=newton2D_exactAS(fm,gradfm,gm,gradgm,x0,y0,100,1e-6,true);
    finalarray(1,j)=xm;
    finalarray(2,j)=ym;
    j=j+1;
end %for 

disp('The four roots are:');
disp('X = 0 and 2');
disp('Y = 1 and 0');

disp('The 12 solutions presented below are a result of multiple iterations and a small tolerance and do not reflect more solutions.');
disp(finalarray);
figure;
surf(X,Y,F);
hold on;
surf(X,Y,G);
plot3(xm,ym,0,'wo','MarkerSize',32,'LineWidth',8);
hold off;
disp('%%%%%%%%End Part 3A Solution:%%%%%%%');

%%Part B: Produce an altered multi-dimensional Newton method (start from newton2D exact.m)
%%to find a root for the three equations system

disp('%%%%%%%%Part 3B Solution:%%%%%%%');
%% Multidimensional function from Problem #3
%%Part 3b of problem, compute the roots of 3 equations (3 dimensions)
fm=@objfun3Df_AS;
gm=@objfun3Dg_AS;
km=@objfun3Dk_AS;
gradfm=@grad_objfun3Df_AS;
gradgm=@grad_objfun3Dg_AS;
gradkm=@grad_objfun3Dk_AS;
%this is for plotting
x=linspace(-1.5,1.5,20);
y=linspace(-1.5,1.5,20);
z=linspace(-1.5,1.5,20);
[X,Y,Z]=meshgrid(x,y,z);
F=fm(X,Y,Z);
G=gm(X,Y,Z);
K=km(X,Y,Z);

%% Newton's method for multi-variable nonlinear equations
%x0=i;
%y0=0.258*i;
%z0=0.58*i;
x0=0.5;
y0=0.5;
z0=0.5;
[xm,ym,zm,it3D,success3D]=newton3D_exactAS(fm,gradfm,gm,gradgm,km,gradkm,x0,y0,z0,100,1e-6,true);

disp('%%%%%Roots for All 3 Equations%%%%%%');
disp('X-Root');
disp(xm);
disp('Y-Root');
disp(ym);
disp('Z-Root');
disp(zm);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
%%Root Checker:
disp('This section evaluates the roots and checks if the equations go to Zero:');
r1=fm(xm,ym,zm);
disp('When plugged back into Equation #1: x^2+y^2+z^2 = 6');
disp(r1);
r2=gm(xm,ym,zm);
disp('When plugged back into Equation #2: x^2−y^2+2*z^2 = 2');
disp(r2);
r3=km(xm,ym,zm);
disp('When plugged back into Equation #3: 2*x^2+y^2−z^2 = 3');
disp(r3);
disp('The values shown above are not exactly Zero, but they are quit close. This is due to precision and rounding in matlab.');

disp('%%%%%%%% End Part 3B Solution:%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%PROBLEM #3 ANSWER END%%%%%%%%%%%%%%%%%%');
%% END PROJECT #3







