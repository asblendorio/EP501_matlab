% Problem 2(a) Write a function that takes in a grid of points describing 
% some independent variable (say xi), and a point to which the data are to be 
% interpolated x′ and finds the index i into the array xi such that: xi ≤ x′ ≤ xi+1.

% Problem 2(b) Use the function from part (a) to construct an additional 
% function that works over a 2D grid x, y. I.e. given two grids xi, yj 
% find the indices i, j such that: xi ≤ x′ ≤ xi+1, yj ≤ y′ ≤ yj+1.

function [ix,xi1,iy,yi1] = bilinear(yg,xg,xprime,yprime)
xint=xprime;
nref1=length(xg);
x=zeros(1,nref1);
x(1)=xg(1);
for i=2:nref1
    x(i) = (xg(i)-xint);
    if x(i) >= 0
        if x(i-1) <= 0 
            xi1 = x(i);
            ix = i-1;
        end %if 
    end %if
end %for 

yint=yprime;
nref2=length(yg);
y=zeros(1,nref2);
y(1)=yg(1);
for i=2:nref2
    y(i) = yg(i)-yint;
    if y(i) >= 0
        if y(i-1) <= 0 
            yi1 = y(i);          
            iy = i-1;
        end %if 
    end %if
end %for 

% disp('The input points are ');
% disp(xint);
% disp(yint);
% disp('X ');
% disp(xg(ix));
% disp(xg(ix+1));
% disp('Y');
% disp(yg(iy));
% disp(yg(iy+1));


%% Plotter for 2A and B
% figure(3);
% title('Problem A and B Indices I and J');
% plot(xint,yint,'*','MarkerSize',5,'MarkerFaceColor','red');
% hold on;
% plot(xg(xi),yg(yi), 'o','MarkerSize',5,'MarkerFaceColor','blue');
% hold on;
% plot(xg(xi+1),yg(yi+1),'o','MarkerSize',5,'MarkerFaceColor','blue');
% hold on;
% plot(xg(xi),yg(yi+1),'o','MarkerSize',5,'MarkerFaceColor','blue');
% hold on;
% plot(xg(xi+1),yg(yi), 'o','MarkerSize',5,'MarkerFaceColor','blue');
% hold off
% 
% 
% xvec=xi1(:);
% yvec=yi1(:);
% f = f2D;
% 
% figure(4);
% imagesc(x,y,f);
% axis xy;
% xlabel('x');
% ylabel('y');
% c=colorbar;
% ylabel(c,'f(x,y)')
% hold on;
% plot(xvec,yvec,'w^','MarkerSize',15,'MarkerFaceColor','white');
% plot(xint,yint,'wo','MarkerSize',20,'MarkerFaceColor','white');
% hold off;
% 


end % function 