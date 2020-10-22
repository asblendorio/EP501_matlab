% Problem 3 B through E 
% This is an implementation of the Newton-Horner Method that will perform
% all tasks laid out in problems 3 B through E. 
% 
function[r,Pnew]= polySolver(a,r0,tol,nmax)
nref=length(a)-1;  
r=zeros(nref,1);    % setting up empty array
Pnew=zeros(nref,1);  % setting up empty array
for i=1:1:nref      % 
    ni=0;
    x=r0;
    err=10;
    while ni<nmax && err>=tol % 
        [v,p]=polynomial(a,x);   % pass in the synthetic polynomial division algorithm to perform factoring 
        [dv,~]=polynomial(p,x);   % computes the derivative of the polynomial using output 
        xn=x-v/dv;
        err=abs(xn-x);              
        ni=ni+1;
        x=xn;
    end %while    
    r(i)=x;                 
    Pnew(i)=ni;
    [~,a]=polynomial(a,x); %performs deflation for each polynomial 
      
end %for

end %function 