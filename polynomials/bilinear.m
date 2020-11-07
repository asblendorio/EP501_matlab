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

function [xi,xi1,yi,yi1] = bilinear(xint,yint,yg,ygi,xg,xgi)
xint=2;
nref1=length(xg);
x=zeros(1,nref1);
for i=1:nref1
    x(i) = (xg(i)-xint);
    if x(i) > 0
        if x(i-1) < 0 
            xi1 = x(i);
            xi = i-1;
        end %if 
    end %if
end %for 

yint=1;
nref2=length(yg);
y=zeros(1,nref2);
for i=1:nref2
    y(i) = yg(i)-yint;
    if y(i) > 0
        if y(i-1) < 0 
            yi1 = y(i);          
            yi = i-1;
        end %if 
    end %if
end %for 

disp(yi1);
disp(yi); 
disp(xi1);
disp(xi);

plot(xint,yint,'*','MarkerSize',5,'MarkerFaceColor','red');
hold on;
plot(xg(xi),yg(yi), 'o','MarkerSize',5,'MarkerFaceColor','blue');
hold on;
plot(xg(xi+1),yg(yi+1),'o','MarkerSize',5,'MarkerFaceColor','blue');
hold on;
plot(xg(xi),yg(yi+1),'o','MarkerSize',5,'MarkerFaceColor','blue');
hold on;
plot(xg(xi+1),yg(yi), 'o','MarkerSize',5,'MarkerFaceColor','blue');
hold off

end % function 