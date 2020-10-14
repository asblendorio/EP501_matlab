%% Author: Alec Sblendorio 
%% Due Date: October 15, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% Project #3
% Problem #1: Finding roots of functions lacking a closed form
% Problem #2: Numerical solution for multiple polynomial roots
% Problem #3: Multivariate root finding
%% Import test data for calibration of software

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
disp('%%%%%%%%Part A Solution:%%%%%%%');

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

disp('%%%%%%%%End Part A Solution:%%%%%%%');

disp('%%%%%%%%Part B Solution:%%%%%%%');
%%In Office Hours, Dr. Z alluded to looping over values form 0 to 20 to
%%calculate the first root of the Bessel function. 
maxit=100;


disp('%%%%%%%%End Part B Solution:%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER END%%%%%%%%%%%%%%%%%%');
%% Problem #2: Numerical solution for multiple polynomial roots
disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
%%Part A: Suppose you have a polynomial of known or given order and need to find all of its roots.
%%Write a block of code or a script that uses the exact Newton method (e.g. the function in the course repository)
%%to find all of the real-valued roots of a polynomial. 
%%Assume that you do not need to identify repeated roots (if any). 
%%Test your code on the polynomial used in the book (Eqn. 3.115).

disp('%%%%%%%%Part A Solution:%%%%%%%');

disp('%%%%%%%%End Part A Solution:%%%%%%%');

%%Part B: Produce an altered version of your code to deal with the fact that there are potentially complex roots
%%to your polynomial. Use this code to find all roots, including complex-valued solutions
%%of the following polynomial (Eqn 3.145 in the book). 


disp('%%%%%%%%Part B Solution:%%%%%%%');

disp('%%%%%%%% End Part B Solution:%%%%%%%');


disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER END%%%%%%%%%%%%%%%%%%');

%% Problem #3: Multivariate root finding
disp('%%%%%%%%%%%%%%%%%%PROBLEM #3 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
%%Part A: Use the multi-dimensional Newton method
%%(given the in course repository: newton2D exact.m) to find all four roots of the system

disp('%%%%%%%%Part A Solution:%%%%%%%');
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
for i=-2.5:0.25:2.5
    x0=i;
    y0=i-0.4;
    [xm,ym,it2D,success2D]=newton2D_exactAS(fm,gradfm,gm,gradgm,x0,y0,100,1e-6,true);
%     finalarray(1,j)=xm;
%     finalarray(2,j)=ym;
    j=j+1;
end %for

disp('Solution');
%disp(finalarray);
disp(xm);
disp(ym);
% disp(xm);
% disp(ym);
% disp(it2D);
% disp(success2D);

figure;
surf(X,Y,F);
hold on;
surf(X,Y,G);
plot3(xm,ym,0,'wo','MarkerSize',32,'LineWidth',8);

disp('%%%%%%%%End Part A Solution:%%%%%%%');

%%Part B: Produce an altered multi-dimensional Newton method (start from newton2D exact.m)
%%to find a root for the three equations system

disp('%%%%%%%%Part B Solution:%%%%%%%');
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


disp('%%%%%%%% End Part B Solution:%%%%%%%');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #3 ANSWER END%%%%%%%%%%%%%%%%%%');
%% END PROJECT #3







