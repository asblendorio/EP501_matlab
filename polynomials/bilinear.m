% Problem 2(a) Write a function that takes in a grid of points describing 
% some independent variable (say xi), and a point to which the data are to be 
% interpolated x′ and finds the index i into the array xi such that: xi ≤ x′ ≤ xi+1.

% Problem 2(b) Use the function from part (a) to construct an additional 
% function that works over a 2D grid x, y. I.e. given two grids xi, yj 
% find the indices i, j such that: xi ≤ x′ ≤ xi+1, yj ≤ y′ ≤ yj+1.

% Problem 2(c) Use your results from parts a and b to create a bilinear 
% interpolation function that takes in a sequence of data points {x′k,yk′ } 
% to which data are being interpolated, a grid xi,yj, and a dataset fij 
% that is defined over this grid and produces bilinearly interpolated values 
% of fk at the points {x′k,yk′ }. Write your program so that the input
% points are simply a flat list and not necessarily a 2D grid of points 

function [yi,y1,yi1,xi,x1,xi1] = bilinear(order,yg,ygi,xg,xgi)
nref=length(yg);
y=ygi(order);
for i=1:1:nref-1
    if (y-yg(i)) < 1
       ynew=i;
       yi=yg(i);
       yi1=yg(i+1);
    end %if
end %for 

nref2=length(xg);
x=xgi(order);
for i=1:1:nref2-1
    if (x-xg(i))<1
        xnew=i;
        xi=xg(i);
        xi1=xg(i+1);
    end %if 
end %for     

disp(yi1);
disp(y);
disp(yi);

disp(xi1);
disp(x);
disp(xi);

end % function 