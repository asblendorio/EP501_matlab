%% Author: Alec Sblendorio 
%% REVISED Due Date: December 7, 2020
%% Subject: EP 501 Numerical Methods for Scientists and Engineers 
%% Project #4 REVISION 
%%Problem 1: This problem concerns least squares and data fitting and requires use of the example dataset from the repository, test lsq.mat, which provides data for variables xi, yi, σyi referenced below.
%%Problem #2: This problem concerns bilinear interpolation methods and requires use of the grid (variables xg,yg) and data samples (f2D) from test interp.mat
% Collaboration with Sho Okayama, Dennis Turano, and Bartosz Kwiecisnki
%% Data Input 
load('test_lsq.mat');
load('test_interp.mat');
%% COMMENTS ON REVISION 
%%The previous submission of Project 4 did not have the correct solutions
%%to problems 2 c and d. I have revised it to properly interpolate the
%%given function and test the results against Matlab's bilinear
%%interpolation function (interp2). Shown in this report are two plots. One
%%for the numerically calculated version and another for Matlab's solver. 
%%Working with Dr. Z to debug the bilinear function, the correct
%%interpolation was performed. The function "bilinear_Rev1.m" and this
%%script "Project4_Rev1.m" have the correct solutions for problem 2. 
%%Dr. Z mentioned it is not necessary to re-submit all the code for the
%%project reconciliation so I trimmed the PDF file down to the relevent
%%information. The .zip file will need the ret of the code to run so I left
%%that in. 

%% Problem #2: This problem concerns bilinear interpolation methods and requires use of the grid (variables xg,yg) and data samples (f2D) from test interp.mat.
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

%% Illustration of bilinear interpolation, single interval of interest

f=f2D;
x1=xg;
y1=yg;
% Manually written
% n = length(f(:,1));
% nref = length(f(:,1));
n=numel(xgi);
nref=numel(ygi);
finterp = zeros(512,512);
for ir1 = 1:1:n-1
    for ir2=1:1:nref-1
        xprime=xgi(ir1);
        yprime=ygi(ir2);
        [ix,xi1,iy,yi1] = bilinear(yg,xg,xprime,yprime);
        x=[xg(ix),xg(ix+1)];
        y=[yg(iy),yg(iy+1)];
        [X,Y] = meshgrid(x,y);
%       fvec=f(ir1:ir1+1,ir2:ir2+1);
        fvec=f(iy:iy+1,ix:ix+1);
        xvec=X(:);
        yvec=Y(:);
        M=[ones(4,1),xvec(:),yvec(:),xvec(:).*yvec(:)];
        [Mmod,order]=Gauss_elim(M,fvec(:));
        avec=backsub(Mmod(order,:));
        finterpmanual(ir1,ir2)=avec(1)+avec(2)*xprime+avec(3)*yprime+avec(4)*xprime*yprime;
        finterp(ir1,ir2)=interp2(X,Y,fvec,xprime,yprime);
        
    end %for
end %for
%% Plotter Numerical Version 
figure(2);
pcolor(finterpmanual);
shading flat;
title('Numerical Bilinear Interpolation');

%% Plotter Matlab version
figure(3);
pcolor(finterp);
shading flat;
title('Matlab Bilinear Interpolation');
hold off;

%% Illustrate cubic spline approximations using Matlab functions
x=linspace(-5,5,15);
y=sin(x);
figure(4);
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