% write a function that performs bilinear interpolation for x 
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
        xi=yg(i);
        xi1=yg(i+1);
    end %if 
end %for     

disp(yi1);
disp(y);
disp(yi);

disp(xi1);
disp(x);
disp(xi);

end % function 