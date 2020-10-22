% Problem 3 B , D, and E 
% This is an implementation of the Newton Method found on pg 194-195.
% that will perform all tasks laid out in problems 3 B through E. 
% This expands on the original synthetic division algorithm (SPDA) 
function[r,Pnew]= polySolver(A,r0,nmax)
nref=length(A)-1;  
r=zeros(nref,1);    % setting up empty array
Pnew=zeros(nref,1);  % setting up empty array
for i=1:1:nref      
    ni=0;
    x=r0;
    while ni<nmax   
        [v,p]=polynomial(A,x);   % pass in the synthetic polynomial division algorithm to perform factoring 
        [dv,~]=polynomial(p,x);   % computes the derivative by applying SPDA  
        xn=x-v/dv;              
        ni=ni+1;
        x=xn;
    end %while    
    r(i)=x;                 
    Pnew(i)=ni;
    [~,A]=polynomial(A,x); %performs deflation for each polynomial
end %for

end %function 