% write a function that performs bilinear interpolation for x 
function [x,y] = bilinear(order,yg,ygi,xg,xgi)
nref=length(yg);
y=ygi(order);
for i=1:nref-1
    if y-yg(i)<1
       yi=i;
       y1=yg(i+1);
       yi1=yg(i);
    end %if
end %for 

nref2=length(xg);
x=xgi(order);
for i=1:nref2-1
    if x-xg(i)<1
        xi=i;
        x1=xg(i+1);
        xi1=xg(i);
    end %if 
end %for     
 
disp(yi);
disp(y1);
disp(yi1);
disp(xi);
disp(x1);
disp(xi1);


end % function 