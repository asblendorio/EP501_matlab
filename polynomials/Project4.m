%% Author: Alec Sblendorio 
%% Due Date: November 6, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% Project #4
%%Problem 1: This problem concerns least squares and data fitting and requires use of the example dataset from the repository, test lsq.mat, which provides data for variables xi, yi, σyi referenced below.
%%Problem #2: This problem concerns bilinear interpolation methods and requires use of the grid (variables xg,yg) and data samples (f2D) from test interp.mat
% Collaboration with Sho Okayama, Dennis Turano, and Bartosz Kwiecisnki
%% Data Input 
load('test_lsq.mat');
load('test_interp.mat');
%% Problem #1: Finding roots of functions lacking a closed form
%%Part A: Write a program that performs a linear least squares fit of a set 
%%of data to a polynomial of arbitrary order n. You may use any functions 
%%in the course repository or that you have written for your homework for 
%%solving the required system of equations; however, you may not use the 
%%built-in MATLAB functions for your solution (only to check the results).

%%Part B: Use your fitting program to fit the test data (yi located at independent variable positions xi) 
%%to a line and a quadratic form. Plot your results and the data on the 
%%same axis so they can be easily compared. Test your results against the
%%built-in Matlab functions polyfit and polyval. Compare the error vectors 
%%and residuals for the two fits.

%%Part C: A rigorous way of deciding between preferred functional forms in fits
%%(from some set of options like linear, quadratic, cubic, quartic, etc.)
%%is to define a goodness-of-fit statistic that quantifies how effective a 
%%particular form is a describing a given data set. This goodness statistic 
%%should balance the need to fit the data (by having more parameters (unknowns)
%%in the fit against the fact that of course one can fit a set of data given
%%enough unknowns (a problem referred to as “overfitting the data”). 
%%The simplest and most commonly used goodness of fit statistic is the
%%reduced Chi-squared statistic. Write a function that evaluates χ2ν for a 
%%polynomial fit of order n.

%%Part D: Use your goodness-of-fit statistic to determine whether is best 
%%to fit these data with a linear, quadratic, or cubic polynomial.
%%Show how you reached your decision.

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
disp('This is for linear fit');
soln1=@lesq;
a=soln1(x,ynoisy,sigmay,1);
disp(a);

disp('This is for quadratic fit');
soln2=@lesq;
b=soln2(x,ynoisy,sigmay,2);
disp(b);

disp('This is for cubic fit');
soln3=@lesq;
c=soln3(x,ynoisy,sigmay,3);
disp(c);

%Comments on Answers:
% Inputting the function lesq.m into this script requires I have at least
% one output argument so it can run. Because of that, there are several
% instances of 'extra' answers throughout the execution of this code. 

disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER END%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%PROBLEM #1 ANSWER PLOTS%%%%%%%%%%%%%%%%%%');
%% Problem #2: This problem concerns bilinear interpolation methods and requires use of the grid (variables xg,yg) and data samples (f2D) from test interp.mat.
%%Part A:  Write a function that takes in a grid of points describing some independent variable
%%(say xi), and a point to which the data are to be interpolated x′ and 
%%finds the index i into the array xi such that: xi ≤ x′ ≤ xi+1.

%%Part B: Use the function from part (a) to construct an additional function 
%%that works over a 2D grid x, y. I.e. given two grids xi, yj find the indices i, j 
%%such that: xi ≤ x′ ≤ xi+1, yj ≤ y′ ≤ yj+1.

%%Part C: Use your results from parts a and b to create a bilinear interpolation function
%%that takes in a sequence of data points {x′k,yk′ } to which data are being 
%%interpolated, a grid xi,yj, and a dataset fij that is defined over this grid
%%and produces bilinearly interpolated values of fk at the points {x′k,yk′ }.
%%Write your program so that the input points are simply a flat list and not
%%necessarily a 2D grid of points (you can always reshape the results later if needed).

%%Part D: Test your results against Matlab’s bilinear interpolation function (interp2)
%%and show that you get the same result. Use the test data from the repository 
%%(test interp.mat). The source grid data are stored in xg,yg, while the value 
%%of the function at those points is in f2D. xgi,ygi are the densely sampled 
%%grid point to which the data are to be interpolated for this test.
disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER BEGIN%%%%%%%%%%%%%%%%%%');
%% Comments on Problem #2A and 2B
% Figure 3 presents xi, yj with the indices i, j indicated with circles and
% asterisks. Or so I believe. I was expecting to inquire about this plot in
% the office hours on November 6, 2020 at 5 PM.
% A point of interest is the values that are output as seen below. The
% bilinear interpolation method was completed for a set of data with a
% hardcoded value (as in xint and yint), and it works. However, when it is
% applied to the test data (test_interp.m) then trouble ensues. 
%% Comments on Problem #2C
% Figure 4 shows what the algorithm produces for a bilinear interpolation function 
% that takes in a sequence of data points {x′k,yk′} to which data are being
% interpolated, a grid xi,yj, and a dataset fij that is defined over this 
% grid and produces bilinearly interpolated values of fk at the points {x′k,yk′}. 
% I do not believe this is the right answer however, this was another
% question I would have asked during the office hours. Trying to debug this
% code while not fully understanding the mathematical implications of it is
% difficult. While I was able to figure out some of the problem, I don't
% quite understand why the code is not interpolating correctly. 
%% Comments on Problem #2D
% Trying to use interp2 I was given an error message as seen below: 
% Error using griddedInterpolant
% The sample points arrays must have the same size as the sample values array.
% 
% Error in interp2>makegriddedinterp (line 228)
%     F = griddedInterpolant(varargin{:});
% 
% Error in interp2 (line 136)
%         F = makegriddedinterp(X, Y, V, method,extrap);
% 
% Error in Project4 (line 126)
% finterp=interp2(X,Y,f,x1,y1);
% 
% I have not yet been able to figure out exactly how the interp2 function
% actually operates so this message is cryptic and discussing with others
% in the class, I was not able to debug it in time. 


%% Illustration of bilinear interpolation, single interval of interest
[xi,xi1,yi,yi1] = bilinear(yg,xg,f2D);
x=[xi,xi1];
y=[yi,yi1];
f=f2D;
x1=xg;
y1=yg;
% Manually written
n = length(f(:,1));
nref = length(f(:,1));

for ir1 = 1:1:n-1
    for ir2=1:1:nref-1
        [X,Y] = meshgrid(x,y);
        fvec=f(ir1:ir1+1,ir2:ir2+1);
        xvec=X(:);
        yvec=Y(:);
        M=[ones(4,1),xvec(:),yvec(:),xvec(:).*yvec(:)];
        [Mmod,order]=Gauss_elim(M,fvec(:));
        avec=backsub(Mmod(order,:));
        finterpmanual(ir1,ir2)=avec(1)+avec(2)*x1(ir1)+avec(3)*y1(ir2)+avec(4)*x1(ir1)*y1(ir2);
    end %for
end %for

% Matlab version
% finterp=interp2(X,Y,f,x1,y1);
% disp('Matlab,GNU/Octave built-in solution:');
% disp(finterp);

%% Illustrate cubic spline approximations using Matlab functions
x=linspace(-5,5,15);
y=sin(x);
figure(2);
plot(x,y,'o','MarkerSize',20);

splinedef=spline(x,y);
x2=linspace(min(x),max(x),256);
y2=ppval(splinedef,x2);
hold on;
plot(x2,y2,'.');

y2true=sin(x2);
plot(x2,y2true,'-');

disp('%%%%%%%%%%%%%%%%%%PROBLEM #2 ANSWER END%%%%%%%%%%%%%%%%%%');

%% END PROJECT #4